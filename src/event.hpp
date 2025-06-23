#ifndef EVENT_HPP
#define EVENT_HPP

#include <random>
#include <vector>
#include <unordered_map>
#include "polyomino.hpp"

#define ATTACHMENT 0
#define DETACHMENT 1

struct StochasticEvent
{
	int type;
	Polyomino polyomino;

	double dt;
	
	StochasticEvent(int type, Polyomino polyomino, double dt)
		: type(type), polyomino(polyomino), dt(dt) {}
};

struct EventBucket
{
	std::vector<Polyomino> polyominoes;
	std::unordered_map<Polyomino, unsigned> indexMap;
	
	virtual double getLogProbability() const = 0;
	virtual int getType() const = 0;

	bool addEvent(Polyomino polyomino);
	bool removeEvent(Polyomino polyomino);
	Polyomino pop(int idx);
};

struct AttachmentEventBucket : public EventBucket
{
	double Gmc, concentration;

	AttachmentEventBucket(double Gmc, double concentration)
		: Gmc(Gmc), concentration(concentration) {}
	
	double getLogProbability() const override;
	int getType() const override;
};

struct DetachmentEventBucket : public EventBucket
{
	double Gse;
	int strength;

	DetachmentEventBucket(double Gse, int strength) : Gse(Gse), strength(strength) {}
	
	double getLogProbability() const override;
	int getType() const override;
};

struct StochasticEventSet
{
	std::vector<EventBucket*> eventBuckets;
	std::vector<AttachmentEventBucket*> attachmentBuckets;
	std::vector<DetachmentEventBucket*> detachmentBuckets;

	StochasticEventSet(const std::vector<PolyominoType*> &types,
			   int maxStrength, double Gmc, double Gse, double kf);
	
	bool addAttachmentEvent(Polyomino polyomino);
	bool addDetachmentEvent(Polyomino polyomino, int strength);
	bool removeAttachmentEvent(Polyomino polyomino);
	bool removeDetachmentEvent(Polyomino polyomino, int strength);

	StochasticEvent sample();
	void print() const;

	double Gmc, Gse, kf;
	
	std::mt19937 gen;
	std::extreme_value_distribution<double> gumbel;
};

#endif
