#include <iostream>
#include <chrono>
#include <cmath>
#include <stdexcept>

#include "stochastic.hpp"

Polyomino EventBucket::pop(int idx)
{
	Polyomino selected = this->polyominoes[idx];
	removeEvent(selected);
	return selected;
}

double AttachmentEventBucket::getLogProbability() const
{
	return log(polyominoes.size() * concentration) - Gmc;
}

double DetachmentEventBucket::getLogProbability() const
{
	return log(polyominoes.size()) - Gse * strength;
}

int AttachmentEventBucket::getType() const { return ATTACHMENT; }
int DetachmentEventBucket::getType() const { return DETACHMENT; }

StochasticEventSet::StochasticEventSet(const std::vector<PolyominoType*> &types,
				       int maxStrength, double Gmc, double Gse, double kf)
	: Gmc(Gmc), Gse(Gse), kf(kf)
{
	for (PolyominoType *type : types) {
		auto *bucket = new AttachmentEventBucket(Gmc, type->concentration);
		eventBuckets.push_back(bucket);
		attachmentBuckets.push_back(bucket);
	}

	for (int i = 0; i <= maxStrength; ++i) {
		auto *bucket = new DetachmentEventBucket(Gse, i);
		eventBuckets.push_back(bucket);
		detachmentBuckets.push_back(bucket);
	}

	this->gen.seed(std::chrono::system_clock::now().time_since_epoch().count());
}

void StochasticEventSet::print() const
{
	std::cout << "stochastic event set===============================" << std::endl;
	for (EventBucket *bucket : eventBuckets) {
		std::cout << "  bucket of type " << bucket->getType() <<
			" with " << bucket->polyominoes.size() << " events: log probability = " <<
			bucket->getLogProbability() << std::endl;
		
		for (Polyomino p : bucket->polyominoes) {
			std::cout << p.type->name << p.offset.toString() << ", ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

StochasticEvent StochasticEventSet::sample()
{
	EventBucket *selected = nullptr;
	double maxValue, g, lp;

	double rateAny = 0.0;

	double logkf = log(kf);
	
	for (EventBucket *bucket : this->eventBuckets) {
		if (!bucket->polyominoes.empty()) {
			lp = bucket->getLogProbability();
			g = lp + gumbel(gen);
			if (lp + logkf > -300)
			{
				rateAny += exp(lp + logkf);
			}
			if (selected == nullptr or g > maxValue) {
				maxValue = g;
				selected = bucket;
			}
		}
	}
	
	if (selected == nullptr)
		throw std::runtime_error("No events to be sampled");
	
	int size = selected->polyominoes.size() - 1;
	auto uniformInt = std::uniform_int_distribution<int>(0, size);
	Polyomino polyomino = selected->polyominoes[uniformInt(gen)];

	auto uniformDouble = std::uniform_real_distribution<double>(0, 1);
	double dt = -log(1-uniformDouble(gen)) / rateAny;
	
	return StochasticEvent(selected->getType(), polyomino, dt);
}

bool EventBucket::addEvent(Polyomino polyomino) {
	if (indexMap.count(polyomino) != 0) {
		// throw std::runtime_error("attemped to add existing event of type "+std::to_string(getType())+": "+polyomino.type->name+" "+polyomino.offset.toString());
		return false;
	}

	indexMap[polyomino] = polyominoes.size();
	polyominoes.push_back(polyomino);
	return true;
}

bool EventBucket::removeEvent(Polyomino polyomino) {
	if (indexMap.count(polyomino) == 0) {
		// throw std::runtime_error("attemped to remove non-existing event of type "+std::to_string(getType())+": "+polyomino.type->name+" "+polyomino.offset.toString());
		return false;
	}

	unsigned index = indexMap[polyomino];
	Polyomino selected = polyominoes[index];
	
	if (index != polyominoes.size() - 1) {
		Polyomino last = polyominoes.back();
		polyominoes[index] = last;
		indexMap[last] = index;
	}

	indexMap.erase(selected);
	polyominoes.pop_back();

	return true;
}

bool StochasticEventSet::addAttachmentEvent(Polyomino polyomino)
{
	return attachmentBuckets[polyomino.type->id]->addEvent(polyomino);
}

bool StochasticEventSet::addDetachmentEvent(Polyomino polyomino, int strength)
{
	return detachmentBuckets[strength]->addEvent(polyomino);
}

bool StochasticEventSet::removeAttachmentEvent(Polyomino polyomino)
{
	return attachmentBuckets[polyomino.type->id]->removeEvent(polyomino);
}

bool StochasticEventSet::removeDetachmentEvent(Polyomino polyomino, int strength)
{
	return detachmentBuckets[strength]->removeEvent(polyomino);
}
