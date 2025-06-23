#include <iostream>
#include <algorithm>
#include <chrono>
#include "stochastic.hpp"

int relStrength(PolyominoType *pt1, PolyominoType *pt2, const Vec3 &offset)
{
    int strength = 0;
    for (const Block &b1 : pt1->blocks) {
	for (const Block &b2 : pt2->blocks) {
	    strength += blockStrength(b1, b2, offset);
	}
    }

    return strength;
}

int PrecomputedMaps::getBinding(PolyominoType *pt1, PolyominoType *pt2, const Vec3 &offset) {
    int id1 = pt1->id;
    int id2 = pt2->id;
    
    std::shared_mutex &mutex = this->bindingMapLocks[id1][id2];

    mutex.lock_shared();
    int strength = -1;
    if (this->bindingMap[id1][id2].count(offset) != 0) {
	strength = this->bindingMap[id1][id2][offset];
    }
    mutex.unlock_shared();

    if (strength >= 0) return strength;

    strength = relStrength(pt1, pt2, offset);

    mutex.lock();
    this->bindingMap[id1][id2][offset] = strength;
    mutex.unlock();

    return strength;
}

bool complementaryDomainLabels(const std::string &a, const std::string &b)
{
    return a + "*" == b or a == b + "*";
}

int blockStrength(const Block &b1, const Block &b2, const Vec3 &offset)
{
    Vec3 diff = (b2.coord + offset) - b1.coord;
    if (std::abs(diff.x) + std::abs(diff.y) + std::abs(diff.z) != 1)
	return 0;

    int strength = 0;
    for (const Domain &d1 : b1.domains) {
	for (const Domain &d2 : b2.domains) {
	    if (d1.strength != d2.strength)
		continue;
	    if (d1.direction + d2.direction != 5)
		continue;
	    if (!complementaryDomainLabels(d1.label, d2.label))
		continue;
	    if (diff == directionOffset(d1.direction))
		strength += d1.strength;
	}
    }

    return strength;
}

AbstractSimulation::AbstractSimulation(
    PrecomputedMaps *sharedData,
    const std::unordered_set<Polyomino> &seedPolyominoes,
    const std::unordered_set<Polyomino> &assemblyPolyominoes,
    SimParams simParams) :
    Simulation(sharedData, seedPolyominoes, assemblyPolyominoes, simParams),
    // int bindingThreshold,
    // DimensionRestriction dimRestrictions, 
    // int outputInterval, std::string outputDir, std::string outputName, int numSteps, int threadId, bool reportNonDeterminism,
    // long long rngSeed)
    bindingThreshold(simParams.bindingThreshold),
    reportNonDeterminism(simParams.reportNonDeterminism)
{
    this->sharedData = sharedData;

    this->steps = 0;

    if (this->rngSeed == -1 or this->threadId != 0) {
	this->rngSeed = std::chrono::system_clock::now().time_since_epoch().count();
    }

    this->gen.seed(this->rngSeed);
    std::cout << "TID " << this->threadId << " using RNG seed: " << this->rngSeed << std::endl;

    for (Polyomino poly : seedPolyominoes)
	this->addPolyomino(poly, true);

    for (Polyomino poly : assemblyPolyominoes)
	this->addPolyomino(poly, false);

    this->initSystemXML();
}

void AbstractSimulation::addFrontier(Polyomino polyomino)
{
    if (!frontier.count(polyomino))
    {
	frontierIds.emplace(polyomino, frontierList.size());
	frontier.insert(polyomino).first;
	// Polyomino frontierPoly = *(frontier.find(polyomino));

	frontierList.push_back(polyomino);

	// add each grid location of the polyomino to our count
	for (Block block : polyomino.type->blocks) {
	    Vec3 coord = block.coord + polyomino.offset;

	    if (this->frontierGridLocationCounts.count(coord) == 0)
		this->frontierGridLocationCounts[coord] = 0;

	    this->frontierGridLocationCounts[coord] += 1;
	}
    }
    else
    {
	std::cout << "Attempting to add existing frontier poly" << std::endl;
    }
}

void AbstractSimulation::removeFrontier(Polyomino polyomino)
{
    if (frontier.count(polyomino))
    {
	int id = frontierIds[polyomino];

	if (id + 1 != frontierList.size())
	{
	    Polyomino curr = polyomino;

	    int lastId = frontierList.size() - 1;
	    Polyomino last = frontierList[lastId];

	    frontierList[lastId] = curr;
	    frontierList[id] = last;

	    frontierIds[last] = id;
	}

	frontier.erase(polyomino);
	frontierIds.erase(polyomino);
	frontierList.pop_back();

	for (Block block : polyomino.type->blocks) {
	    Vec3 coord = block.coord + polyomino.offset;

	    int &count = this->frontierGridLocationCounts[coord];
	    count -= 1;
	    if (count <= 0) {
		this->frontierGridLocationCounts.erase(coord);
	    }
	}
    }
    else
    {
	std::cout << "Attempting to remove non-existing frontier poly" << std::endl;
    }
}

void AbstractSimulation::checkNonDeterminism() {
    for (auto it = this->frontierGridLocationCounts.begin(); it != this->frontierGridLocationCounts.end(); ++it) {
	Vec3 coord = it->first;
	int count = it->second;

	// non-determinism found
	if (count > 1) {
	    std::cout << "TID " << this->threadId << ": Non-Determinism detected at coordinates " << coord.toString() << std::endl;
	}
    }
}

void AbstractSimulation::addPolyomino(Polyomino polyomino, bool seed)
{
    int id1 = polyomino.type->id;
	
    // update neighbors
    ShapeType *sType = this->sharedData->polyominoShapes[id1];
    Shape shape(sType, polyomino.offset);
	
    if (locationMap.count(shape)) {
	std::string errString = "attempting to add polyomino \"" +
	    polyomino.type->name + "\" in occupied location " +
	    polyomino.offset.toString();
		
	throw std::runtime_error(errString);
    }

    int strength = strengthMap[polyomino];

    polyomino.stepAdded = this->steps;

    polyomino.initialBindStrength = strength;
    polyomino.currentBindStrength = strength;

    locationMap.emplace(shape, polyomino);
    if (seed) seedSet.insert(polyomino);
    // else eventSet->addDetachmentEvent(polyomino, strength);

    for (auto pair : this->sharedData->potentialBindingList[id1]) {
	PolyominoType *pType = pair.first;
	ShapeType *nbrType = this->sharedData->polyominoShapes[pType->id];

	Vec3 offset = pair.second;
	Shape nbr(nbrType, offset + polyomino.offset);

	// ignore neighbors that don't satisfy dim restrictions
	if (nbr.offset < dimRestriction.min or nbr.offset >= dimRestriction.max)
	    continue;

	// polyomino already exists here, skip
	if (locationMap.count(nbr))
	    continue;

	Polyomino polyominoNbr(pType, nbr.offset);
	int id2 = pType->id;

	int relStrength = this->sharedData->getBinding(polyomino.type, pType, offset);
	
	int &strengthNbr = strengthMap[polyominoNbr];
	int oldStrengthNbr = strengthNbr;
	strengthNbr += relStrength;

	if (overlappedMap[nbr] == 0 and
	    oldStrengthNbr < this->bindingThreshold and 
	    strengthNbr >= this->bindingThreshold) {
	    addFrontier(polyominoNbr);
	}
    }

    /*
    for (Shape nbr : this->sharedData->neighborLists[sType->id]) {

	Vec3 offset = nbr.offset;
	nbr.offset += polyomino.offset;

	// ignore neighbors that don't satisfy dim restrictions
	if (nbr.offset < dimRestriction.min or nbr.offset >= dimRestriction.max)
	    continue;

	// polyomino already exists here, skip
	if (locationMap.count(nbr))
	    continue;

	for (PolyominoType *pType : nbr.type->polyominoTypes) {
	    Polyomino polyominoNbr(pType, nbr.offset);

	    int id2 = pType->id;
	    int relStrength = this->sharedData->getBinding(polyomino.type, pType ,offset);
            
	    int &strengthNbr = strengthMap[polyominoNbr];
	    int oldStrengthNbr = strengthNbr;
	    strengthNbr += relStrength;

	    if (overlappedMap[nbr] == 0 and
		oldStrengthNbr < this->bindingThreshold and 
		strengthNbr >= this->bindingThreshold) {
		addFrontier(polyominoNbr);
	    }
	}
    }
    */

    // remove overlapping attachment events
    for (Shape ovr : this->sharedData->overlapLists[sType->id]) {
	ovr.offset += polyomino.offset;

	if (ovr.offset < dimRestriction.min or ovr.offset >= dimRestriction.max)
	    continue;
		
	int overlappedOvr = ++(overlappedMap[ovr]);

	for (PolyominoType *pType : ovr.type->polyominoTypes) {
	    Polyomino polyominoOvr(pType, ovr.offset);
	    int strengthOvr = strengthMap[polyominoOvr];
	    if (strengthOvr >= this->bindingThreshold and overlappedOvr == 1) {
		// frontier.erase(polyominoOvr);
		removeFrontier(polyominoOvr);
	    }
	}
    }
}

bool AbstractSimulation::step()
{
    if (this->reportNonDeterminism)
	this->checkNonDeterminism();

    int size = this->frontierList.size();
    if (size == 0)
	return false;

    auto uniformInt = std::uniform_int_distribution<int>(0, size-1);

    Polyomino polyomino = this->frontierList[uniformInt(gen)];
    // this->removeFrontier(polyomino);
    this->addPolyomino(polyomino);
    this->steps += 1;

    return true;
}

void AbstractSimulation::operator()()
{
    for (int i = 1; i <= this->numSteps; i++)
    {
	if(!this->step())
	{
	    std::cout << "THREAD " << this->threadId << ": ";
	    std::cout << "empty frontier, step " << i << "/" << this->numSteps << " (" << this->locationMap.size() << " tiles)" << std::endl;
	    this->writeSystemXML();
	    break;
	}
	if (i % this->outputInterval == 0)
	{
	    std::cout << "THREAD " << this->threadId << ": ";
	    std::cout << "step " << i << "/" << this->numSteps << " (" << this->locationMap.size() << " tiles)" << std::endl;
	    this->writeSystemXML();
	}
    }
}

void AbstractSimulation::initSystemXML()
{
    this->systemXMLDocument = new pugi::xml_document();

    pugi::xml_node sysNode = this->systemXMLDocument->append_child("PolyominoSystem");

    // todo
    // pugi::xml_node dimRestrictionsNode = sysNode.append_child("dim_restrictions");
    // dimRestrictionsNode.append_child(pugi::node_pcdata).set_value("");

    pugi::xml_node rotationNode = sysNode.append_child("rotation");
    rotationNode.append_child(pugi::node_pcdata).set_value("False");

    pugi::xml_node minBindingNode = sysNode.append_child("bindingThreshold");
    minBindingNode.append_child(pugi::node_pcdata).set_value(std::to_string(this->bindingThreshold).c_str());
	
    pugi::xml_node stepsNode = sysNode.append_child("steps");
    stepsNode.append_child(pugi::node_pcdata).set_value(std::to_string(this->steps).c_str());

    pugi::xml_node typesNode = sysNode.append_child("PolyominoTypes");
    for (const PolyominoType *type : this->sharedData->polyominoTypes)
    {
	pugi::xml_node polyominoNode = typesNode.append_child("PolyominoType");
		
	pugi::xml_node nameNode = polyominoNode.append_child("name");
	nameNode.append_child(pugi::node_pcdata).set_value(type->name.c_str());

	pugi::xml_node colorNode = polyominoNode.append_child("color");
	colorNode.append_child(pugi::node_pcdata).set_value(type->color.c_str());

	pugi::xml_node concNode = polyominoNode.append_child("concentration");
	concNode.append_child(pugi::node_pcdata).set_value(std::to_string(type->concentration).c_str());

	pugi::xml_node blocksNode = polyominoNode.append_child("blocks");
	for (const Block &block : type->blocks)
	{
	    pugi::xml_node blockNode = blocksNode.append_child("block");

	    pugi::xml_node coordsNode = blockNode.append_child("coords");
	    coordsNode.append_child(pugi::node_pcdata).set_value(block.coord.toString().c_str());

	    pugi::xml_node domainsNode = blockNode.append_child("domains");
	    for (const Domain &domain : block.domains)
	    {
		pugi::xml_node domainNode = domainsNode.append_child("domain");

		pugi::xml_node labelNode = domainNode.append_child("label");
		labelNode.append_child(pugi::node_pcdata).set_value(domain.label.c_str());

		pugi::xml_node dirNode = domainNode.append_child("direction");
		char dirStr[2] = " ";
		switch(domain.direction)
		{
		case NORTH: dirStr[0] = 'N'; break;
		case SOUTH: dirStr[0] = 'S'; break;
		case EAST:  dirStr[0] = 'E'; break;
		case WEST:  dirStr[0] = 'W'; break;
		case UP:    dirStr[0] = 'U'; break;
		case DOWN:  dirStr[0] = 'D'; break;
		}
		dirNode.append_child(pugi::node_pcdata).set_value(dirStr);

		pugi::xml_node strengthNode = domainNode.append_child("strength");
		strengthNode.append_child(pugi::node_pcdata).set_value(std::to_string(domain.strength).c_str());
	    }
	}
    }

    pugi::xml_node seedNode = sysNode.append_child("seed");

    pugi::xml_node seedPolysNode = seedNode.append_child("Polyominoes");
    for (const Polyomino &poly : this->seedSet)
    {
	pugi::xml_node polyNode = seedPolysNode.append_child("Polyomino");
	pugi::xml_node typeNode = polyNode.append_child("PolyominoType");
	typeNode.append_child(pugi::node_pcdata).set_value(poly.type->name.c_str());

	pugi::xml_node transNode = polyNode.append_child("translation");
	transNode.append_child(pugi::node_pcdata).set_value(poly.offset.toString().c_str());
    }

    pugi::xml_node assemblyNode = sysNode.append_child("assembly");

    pugi::xml_node assemblyPolysNode = assemblyNode.append_child("Polyominoes");
    for (const auto &kv : this->locationMap)
    {
	const Polyomino &poly = kv.second;
	pugi::xml_node polyNode = assemblyPolysNode.append_child("Polyomino");
	pugi::xml_node typeNode = polyNode.append_child("PolyominoType");
	typeNode.append_child(pugi::node_pcdata).set_value(poly.type->name.c_str());

	pugi::xml_node transNode = polyNode.append_child("translation");
	transNode.append_child(pugi::node_pcdata).set_value(poly.offset.toString().c_str());
    }


    // sysNode.remove_child(systemXMLDocument.child("assembly"));
    // pugi::xml_node assemblyNode = sysNode.append_child("assembly");
    // pugi::xml_node polyominoesNode = assemblyNode.append_child("Polyominoes");

    // for (auto kv : sim->locationMap) {
    // 	Polyomino p = kv.second;

    // 	pugi::xml_node polyominoNode = polyominoesNode.append_child("Polyomino");
    // 	pugi::xml_node typeNode = polyominoNode.append_child("PolyominoType");
    // 	typeNode.append_child(pugi::node_pcdata).set_value(p.type->name.c_str());
    // 	pugi::xml_node translationNode = polyominoNode.append_child("translation");
    // 	translationNode.append_child(pugi::node_pcdata).set_value(p.offset.toString().c_str());
    // }
}

void AbstractSimulation::writeSystemXML(std::string note)
{
    pugi::xml_node sysNode = this->systemXMLDocument->child("PolyominoSystem");

    sysNode.remove_child("steps");
    pugi::xml_node stepsNode = sysNode.prepend_child("steps");
    stepsNode.append_child(pugi::node_pcdata).set_value(std::to_string(this->steps).c_str());

    sysNode.remove_child("assembly");
    pugi::xml_node assemblyNode = sysNode.append_child("assembly");

    pugi::xml_node assemblyPolysNode = assemblyNode.append_child("Polyominoes");
    for (const auto &kv : this->locationMap)
    {
	const Polyomino &poly = kv.second;
	pugi::xml_node polyNode = assemblyPolysNode.append_child("Polyomino");
	pugi::xml_node typeNode = polyNode.append_child("PolyominoType");
	typeNode.append_child(pugi::node_pcdata).set_value(poly.type->name.c_str());

	pugi::xml_node transNode = polyNode.append_child("translation");
	transNode.append_child(pugi::node_pcdata).set_value(poly.offset.toString().c_str());

	pugi::xml_node polyStepNode = polyNode.append_child("stepAdded");
	polyStepNode.append_child(pugi::node_pcdata).set_value(std::to_string(poly.stepAdded).c_str());

	// pugi::xml_node polyTimeNode = polyNode.append_child("timeAdded");
	// polyTimeNode.append_child(pugi::node_pcdata).set_value(std::to_string(poly.timeAdded).c_str());

	pugi::xml_node polyBindNode = polyNode.append_child("bindStrength");
	polyBindNode.append_child(pugi::node_pcdata).set_value(std::to_string(poly.currentBindStrength).c_str());

	// pugi::xml_node polyInitialBindNode = polyNode.append_child("initialBindStrength");
	// polyInitialBindNode.append_child(pugi::node_pcdata).set_value(std::to_string(poly.initialBindStrength).c_str());
    }
    // for (const auto &kv : this->locationMap)
    // {
    // 	const Polyomino &poly = kv.second;
    // 	pugi::xml_node polyNode = assemblyPolysNode.append_child("Polyomino");
    // 	pugi::xml_node typeNode = polyNode.append_child("PolyominoType");
    // 	typeNode.append_child(pugi::node_pcdata).set_value(poly.type->name.c_str());

    // 	pugi::xml_node transNode = polyNode.append_child("translation");
    // 	transNode.append_child(pugi::node_pcdata).set_value(poly.offset.toString().c_str());
    // }

    std::string filename = this->outputDir + this->outputName;
    filename += "_tid" + std::to_string(this->threadId);
    filename += "_step" + std::to_string(this->steps);
    filename += ".xml";

    this->systemXMLDocument->save_file(filename.c_str());

}

KineticSimulation::KineticSimulation(PrecomputedMaps *sharedData,
				     const std::unordered_set<Polyomino> &seedPolyominoes,
				     const std::unordered_set<Polyomino> &assemblyPolyominoes,
				     SimParams simParams) :
    // double Gmc, double Gse, double kf,
    // DimensionRestriction dimRestrictions, int minBinding, 
    // int outputInterval, std::string outputDir, std::string outputName, int numSteps, int maxAttachments, int threadId, long long rngSeed)
    Simulation(sharedData, seedPolyominoes, assemblyPolyominoes, simParams),
    Gmc(simParams.Gmc), Gse(simParams.Gse), kf(simParams.kf),
    minBinding(simParams.minBinding),
    maxAttachments(simParams.maxAttachments)
{
    this->sharedData = sharedData;

    int maxStrength = 0;
    for (PolyominoType *type : sharedData->polyominoTypes)
    {
	// std::cout << type->name << " " << std::endl;
	int strength = 0;
	for (const Block &block : type->blocks)
	    for (const Domain &domain : block.domains)
		strength += domain.strength;

	if (maxStrength < strength)
	    maxStrength = strength;
    }

    if (this->rngSeed == -1 or this->threadId != 0) {
	this->rngSeed = std::chrono::system_clock::now().time_since_epoch().count();
    }

    this->eventSet = new StochasticEventSet(sharedData->polyominoTypes, maxStrength, Gmc, Gse, kf);
    this->eventSet->gen.seed(this->rngSeed);

    std::cout << "TID " << this->threadId << " using RNG seed: " << this->rngSeed << std::endl;

    this->time = 0.0;
    this->steps = 0;

    for (Polyomino poly : seedPolyominoes)
	this->addPolyomino(poly, true);

    for (Polyomino poly : assemblyPolyominoes)
	this->addPolyomino(poly, false);

    this->initSystemXML();
}

void KineticSimulation::step()
{
    StochasticEvent event = eventSet->sample();
    if (event.type == ATTACHMENT) {
	addPolyomino(event.polyomino);
    } else {
	removePolyomino(event.polyomino);
    }
    std::cout << "C" << std::endl;
    this->time += event.dt;
    this->steps += 1;
}

void KineticSimulation::operator()()
{
    std::string note = "";
    int i;

    for (i = 1; i <= this->numSteps; i++)
    {
	this->step();
	
	if (this->maxAttachments > 0) {
	    if (this->locationMap.size() - this->seedSet.size() >= this->maxAttachments) {
		std::cout << "THREAD " << this->threadId << ": ";
		std::cout << "MAX ATTACHMENTS REACHED: step " << i << "/" << this->numSteps << " (" << this->locationMap.size() << " tiles)" << std::endl;
		note = "MAX ATTACHMENTS REACHED";
		break;
	    }
	}

	if (i % this->outputInterval == 0)
	{
	    std::cout << "THREAD " << this->threadId << ": ";
	    std::cout << "step " << i << "/" << this->numSteps << " (" << this->locationMap.size() << " tiles)" << std::endl;
			
	    if (i != this->numSteps) {
		this->writeSystemXML();
	    }
	}
    }

    this->writeSystemXML(note);
}

void KineticSimulation::initSystemXML()
{
    this->systemXMLDocument = new pugi::xml_document();

    pugi::xml_node sysNode = this->systemXMLDocument->append_child("PolyominoSystem");

    pugi::xml_node GmcNode = sysNode.append_child("Gmc");
    GmcNode.append_child(pugi::node_pcdata).set_value(std::to_string(this->Gmc).c_str());

    pugi::xml_node GseNode = sysNode.append_child("Gse");
    GseNode.append_child(pugi::node_pcdata).set_value(std::to_string(this->Gse).c_str());

    pugi::xml_node kfNode = sysNode.append_child("forwardRate");
    kfNode.append_child(pugi::node_pcdata).set_value(std::to_string(this->kf).c_str());

    // todo
    // pugi::xml_node dimRestrictionsNode = sysNode.append_child("dim_restrictions");
    // dimRestrictionsNode.append_child(pugi::node_pcdata).set_value("");

    pugi::xml_node rotationNode = sysNode.append_child("rotation");
    rotationNode.append_child(pugi::node_pcdata).set_value("False");

    pugi::xml_node minBindingNode = sysNode.append_child("min_binding");
    minBindingNode.append_child(pugi::node_pcdata).set_value(std::to_string(this->minBinding).c_str());

    pugi::xml_node maxAttachmentsNode = sysNode.append_child("max_attachments");
    maxAttachmentsNode.append_child(pugi::node_pcdata).set_value(std::to_string(this->maxAttachments).c_str());

    pugi::xml_node timeNode = sysNode.append_child("time");
    timeNode.append_child(pugi::node_pcdata).set_value(std::to_string(this->time).c_str());
    // timeNode.append_child(pugi::node_pcdata).set_value(this->time);
	
    pugi::xml_node stepsNode = sysNode.append_child("steps");
    stepsNode.append_child(pugi::node_pcdata).set_value(std::to_string(this->steps).c_str());

    pugi::xml_node typesNode = sysNode.append_child("PolyominoTypes");
    for (const PolyominoType *type : this->sharedData->polyominoTypes)
    {
	pugi::xml_node polyominoNode = typesNode.append_child("PolyominoType");
		
	pugi::xml_node nameNode = polyominoNode.append_child("name");
	nameNode.append_child(pugi::node_pcdata).set_value(type->name.c_str());

	pugi::xml_node colorNode = polyominoNode.append_child("color");
	colorNode.append_child(pugi::node_pcdata).set_value(type->color.c_str());

	pugi::xml_node concNode = polyominoNode.append_child("concentration");
	concNode.append_child(pugi::node_pcdata).set_value(std::to_string(type->concentration).c_str());

	pugi::xml_node blocksNode = polyominoNode.append_child("blocks");
	for (const Block &block : type->blocks)
	{
	    pugi::xml_node blockNode = blocksNode.append_child("block");

	    pugi::xml_node coordsNode = blockNode.append_child("coords");
	    coordsNode.append_child(pugi::node_pcdata).set_value(block.coord.toString().c_str());

	    pugi::xml_node domainsNode = blockNode.append_child("domains");
	    for (const Domain &domain : block.domains)
	    {
		pugi::xml_node domainNode = domainsNode.append_child("domain");

		pugi::xml_node labelNode = domainNode.append_child("label");
		labelNode.append_child(pugi::node_pcdata).set_value(domain.label.c_str());

		pugi::xml_node dirNode = domainNode.append_child("direction");
		char dirStr[2] = " ";
		switch(domain.direction)
		{
		case NORTH: dirStr[0] = 'N'; break;
		case SOUTH: dirStr[0] = 'S'; break;
		case EAST:  dirStr[0] = 'E'; break;
		case WEST:  dirStr[0] = 'W'; break;
		case UP:    dirStr[0] = 'U'; break;
		case DOWN:  dirStr[0] = 'D'; break;
		}
		dirNode.append_child(pugi::node_pcdata).set_value(dirStr);

		pugi::xml_node strengthNode = domainNode.append_child("strength");
		strengthNode.append_child(pugi::node_pcdata).set_value(std::to_string(domain.strength).c_str());
	    }
	}
    }

    pugi::xml_node seedNode = sysNode.append_child("seed");

    pugi::xml_node seedPolysNode = seedNode.append_child("Polyominoes");
    for (const Polyomino &poly : this->seedSet)
    {
	pugi::xml_node polyNode = seedPolysNode.append_child("Polyomino");
	pugi::xml_node typeNode = polyNode.append_child("PolyominoType");
	typeNode.append_child(pugi::node_pcdata).set_value(poly.type->name.c_str());

	pugi::xml_node transNode = polyNode.append_child("translation");
	transNode.append_child(pugi::node_pcdata).set_value(poly.offset.toString().c_str());
    }

    pugi::xml_node assemblyNode = sysNode.append_child("assembly");

    pugi::xml_node assemblyPolysNode = assemblyNode.append_child("Polyominoes");
    for (const auto &kv : this->locationMap)
    {
	const Polyomino &poly = kv.second;
	pugi::xml_node polyNode = assemblyPolysNode.append_child("Polyomino");
	pugi::xml_node typeNode = polyNode.append_child("PolyominoType");
	typeNode.append_child(pugi::node_pcdata).set_value(poly.type->name.c_str());

	pugi::xml_node transNode = polyNode.append_child("translation");
	transNode.append_child(pugi::node_pcdata).set_value(poly.offset.toString().c_str());

	pugi::xml_node polyStepNode = polyNode.append_child("stepAdded");
	polyStepNode.append_child(pugi::node_pcdata).set_value(std::to_string(poly.stepAdded).c_str());

	pugi::xml_node polyTimeNode = polyNode.append_child("timeAdded");
	polyTimeNode.append_child(pugi::node_pcdata).set_value(std::to_string(poly.timeAdded).c_str());

	pugi::xml_node polyBindNode = polyNode.append_child("bindStrength");
	polyBindNode.append_child(pugi::node_pcdata).set_value(std::to_string(poly.currentBindStrength).c_str());

	pugi::xml_node polyInitialBindNode = polyNode.append_child("initialBindStrength");
	polyInitialBindNode.append_child(pugi::node_pcdata).set_value(std::to_string(poly.initialBindStrength).c_str());
    }
}

void KineticSimulation::writeSystemXML(std::string note)
{
    pugi::xml_node sysNode = this->systemXMLDocument->child("PolyominoSystem");

    sysNode.remove_child("time");
    pugi::xml_node timeNode = sysNode.prepend_child("time");
    timeNode.append_child(pugi::node_pcdata).set_value(std::to_string(this->time).c_str());
    // timeNode.append_child(pugi::node_pcdata).set_value(this->time);
	
    sysNode.remove_child("steps");
    pugi::xml_node stepsNode = sysNode.prepend_child("steps");
    stepsNode.append_child(pugi::node_pcdata).set_value(std::to_string(this->steps).c_str());

    sysNode.remove_child("assembly");
    pugi::xml_node assemblyNode = sysNode.append_child("assembly");

    sysNode.remove_child("note");
    pugi::xml_node noteNode = sysNode.prepend_child("note");
    noteNode.append_child(pugi::node_pcdata).set_value(note.c_str());

    pugi::xml_node assemblyPolysNode = assemblyNode.append_child("Polyominoes");
    for (const auto &kv : this->locationMap)
    {
	const Polyomino &poly = kv.second;
	pugi::xml_node polyNode = assemblyPolysNode.append_child("Polyomino");
	pugi::xml_node typeNode = polyNode.append_child("PolyominoType");
	typeNode.append_child(pugi::node_pcdata).set_value(poly.type->name.c_str());

	pugi::xml_node transNode = polyNode.append_child("translation");
	transNode.append_child(pugi::node_pcdata).set_value(poly.offset.toString().c_str());

	pugi::xml_node polyStepNode = polyNode.append_child("stepAdded");
	polyStepNode.append_child(pugi::node_pcdata).set_value(std::to_string(poly.stepAdded).c_str());

	pugi::xml_node polyTimeNode = polyNode.append_child("timeAdded");
	polyTimeNode.append_child(pugi::node_pcdata).set_value(std::to_string(poly.timeAdded).c_str());

	pugi::xml_node polyBindNode = polyNode.append_child("bindStrength");
	polyBindNode.append_child(pugi::node_pcdata).set_value(std::to_string(poly.currentBindStrength).c_str());

	pugi::xml_node polyInitialBindNode = polyNode.append_child("initialBindStrength");
	polyInitialBindNode.append_child(pugi::node_pcdata).set_value(std::to_string(poly.initialBindStrength).c_str());
    }

    std::string filename = this->outputDir + this->outputName;
    filename += "_tid" + std::to_string(this->threadId);
    filename += "_step" + std::to_string(this->steps);
    filename += ".xml";

    this->systemXMLDocument->save_file(filename.c_str());

}

void KineticSimulation::addPolyomino(Polyomino polyomino, bool seed)
{
    int id1 = polyomino.type->id;

    // update neighbors
    ShapeType *sType = this->sharedData->polyominoShapes[id1];
    Shape shape(sType, polyomino.offset);

    if (locationMap.count(shape)) {
	std::string errString = "attempting to add polyomino \"" +
	    polyomino.type->name + "\" in occupied location " +
	    polyomino.offset.toString();
		
	throw std::runtime_error(errString);
    }

    int strength = strengthMap[polyomino];

    polyomino.stepAdded = this->steps;
    polyomino.timeAdded = this->time;

    polyomino.initialBindStrength = strength;
    polyomino.currentBindStrength = strength;

    locationMap.emplace(shape, polyomino);

    if (seed) seedSet.insert(polyomino);
    else eventSet->addDetachmentEvent(polyomino, strength);
    
    if (this->minBinding > 0) {
	for (auto pair : this->sharedData->potentialBindingList[id1]) {
	    PolyominoType *pType = pair.first;
	    ShapeType *nbrType = this->sharedData->polyominoShapes[pType->id];

	    Vec3 offset = pair.second;
	    Shape nbr(nbrType, offset + polyomino.offset);

	    // ignore neighbors that don't satisfy dim restrictions
	    if (nbr.offset < dimRestriction.min or nbr.offset >= dimRestriction.max)
		continue;
	    
	    Polyomino existing = polyomino;
	    if (locationMap.count(nbr))
		existing = locationMap[nbr];
	    
	    Polyomino polyominoNbr(pType, nbr.offset);
	    int id2 = pType->id;

	    int relStrength = this->sharedData->getBinding(polyomino.type, pType, offset);
		
	    int &strengthNbr = strengthMap[polyominoNbr];
	    int oldStrengthNbr = strengthNbr;
	    strengthNbr += relStrength;

	    if (polyominoNbr == existing and seedSet.count(polyominoNbr) == 0) {
		locationMap[nbr].currentBindStrength = strengthNbr;

		eventSet->removeDetachmentEvent(polyominoNbr, oldStrengthNbr);
		eventSet->addDetachmentEvent(polyominoNbr, strengthNbr);
	    } else {
		if (overlappedMap[nbr] == 0 and
		    (minBinding == 0 or oldStrengthNbr < minBinding) and
		    strengthNbr >= minBinding) {
		    eventSet->addAttachmentEvent(polyominoNbr);
		}
	    }
	}
    } else {
	for (Shape nbr : this->sharedData->neighborLists[sType->id]) {
	    Vec3 offset = nbr.offset;
	    nbr.offset += polyomino.offset;

	    // ignore neighbors that don't satisfy dim restrictions
	    if (nbr.offset < dimRestriction.min or nbr.offset >= dimRestriction.max)
		continue;

	    Polyomino existing = polyomino;
	    if (locationMap.count(nbr))
		existing = locationMap[nbr];
		
	    for (PolyominoType *pType : nbr.type->polyominoTypes) {

		Polyomino polyominoNbr(pType, nbr.offset);

		int id2 = pType->id;
		//int relStrength = this->sharedData->bindingMap[id1][id2][offset];
		int relStrength = this->sharedData->getBinding(polyomino.type, pType, offset);
		
		int &strengthNbr = strengthMap[polyominoNbr];
		int oldStrengthNbr = strengthNbr;
		strengthNbr += relStrength;

		if (polyominoNbr == existing and seedSet.count(polyominoNbr) == 0) {
		    locationMap[nbr].currentBindStrength = strengthNbr;

		    eventSet->removeDetachmentEvent(polyominoNbr, oldStrengthNbr);
		    eventSet->addDetachmentEvent(polyominoNbr, strengthNbr);
		} else {
		    if (overlappedMap[nbr] == 0 and
			(minBinding == 0 or oldStrengthNbr < minBinding) and
			strengthNbr >= minBinding) {
			eventSet->addAttachmentEvent(polyominoNbr);
		    }
		}
	    }
	}
    }

    // remove overlapping attachment events
    for (Shape ovr : this->sharedData->overlapLists[sType->id]) {
	ovr.offset += polyomino.offset;

	if (ovr.offset < dimRestriction.min or ovr.offset >= dimRestriction.max)
	    continue;
		
	int overlappedOvr = ++(overlappedMap[ovr]);

	for (PolyominoType *pType : ovr.type->polyominoTypes) {
	    Polyomino polyominoOvr(pType, ovr.offset);
	    int strengthOvr = strengthMap[polyominoOvr];
	    if (strengthOvr >= minBinding and overlappedOvr == 1) {
		// std::cout << polyominoOvr.type->name << " " << polyominoOvr.offset.toString() << std::endl;
		eventSet->removeAttachmentEvent(polyominoOvr);
	    }
	}
    }
}

void KineticSimulation::removePolyomino(Polyomino polyomino)
{
    int id1 = polyomino.type->id;
    // remove state from map and update neighbors
    ShapeType *sType = this->sharedData->polyominoShapes[id1];
    Shape shape(sType, polyomino.offset);
	
    if (!locationMap.count(shape)) {
	std::string errString = "attempting to remove polyomino \"" +
	    polyomino.type->name + "\" from unoccupied location " +
	    polyomino.offset.toString();
	throw std::runtime_error(errString);
    }

    locationMap.erase(shape);
    int strength = strengthMap[polyomino];
    eventSet->removeDetachmentEvent(polyomino, strength);

    if (this->minBinding > 0) {
	for (auto pair : this->sharedData->potentialBindingList[id1]) {
	    PolyominoType *pType = pair.first;
	    ShapeType *nbrType = this->sharedData->polyominoShapes[pType->id];

	    Vec3 offset = pair.second;
	    Shape nbr(nbrType, offset + polyomino.offset);

	    // ignore neighbors that don't satisfy dim restrictions
	    if (nbr.offset < dimRestriction.min or nbr.offset >= dimRestriction.max)
		continue;
	    
	    Polyomino existing = polyomino;
	    if (locationMap.count(nbr))
		existing = locationMap[nbr];
	    
	    Polyomino polyominoNbr(pType, nbr.offset);
	    int id2 = pType->id;

	    int relStrength = this->sharedData->getBinding(polyomino.type, pType, offset);

	    int &strengthNbr = strengthMap[polyominoNbr];
	    int oldStrengthNbr = strengthNbr;
	    strengthNbr -= relStrength;

	    if (polyominoNbr == existing and seedSet.count(polyominoNbr) == 0) {
		eventSet->removeDetachmentEvent(polyominoNbr, oldStrengthNbr);
		eventSet->addDetachmentEvent(polyominoNbr, strengthNbr);
	    } else {
		if (overlappedMap[nbr] == 0 and
		    oldStrengthNbr >= minBinding and
		    strengthNbr < minBinding) {
		    eventSet->removeAttachmentEvent(polyominoNbr);
		}
	    }
	}
    } else {
	for (Shape nbr : this->sharedData->neighborLists[sType->id]) {
	    Vec3 offset = nbr.offset;
	    nbr.offset += polyomino.offset;

	    // ignore neighbors that don't satisfy dim restrictions
	    if (nbr.offset < dimRestriction.min or nbr.offset >= dimRestriction.max)
		continue;

	    Polyomino existing = polyomino;
	    if (locationMap.count(nbr))
		existing = locationMap[nbr];
		
	    for (PolyominoType *pType : nbr.type->polyominoTypes) {
		Polyomino polyominoNbr(pType, nbr.offset);

		int id2 = pType->id;
		//int relStrength = this->sharedData->bindingMap[id1][id2][offset];
		int relStrength = this->sharedData->getBinding(polyomino.type, pType, offset);
	    

		int &strengthNbr = strengthMap[polyominoNbr];
		int oldStrengthNbr = strengthNbr;
		strengthNbr -= relStrength;

		if (polyominoNbr == existing and seedSet.count(polyominoNbr) == 0) {
		    eventSet->removeDetachmentEvent(polyominoNbr, oldStrengthNbr);
		    eventSet->addDetachmentEvent(polyominoNbr, strengthNbr);
		} else {
		    if (overlappedMap[nbr] == 0 and
			oldStrengthNbr >= minBinding and
			strengthNbr < minBinding) {
			eventSet->removeAttachmentEvent(polyominoNbr);
		    }
		}
	    }
	}
    }

    for (Shape ovr : this->sharedData->overlapLists[sType->id]) {
	ovr.offset += polyomino.offset;

	if (ovr.offset < dimRestriction.min or ovr.offset >= dimRestriction.max)
	    continue;
		
	int overlappedOvr = --(overlappedMap[ovr]);

	for (PolyominoType *pType : ovr.type->polyominoTypes) {
	    Polyomino polyominoOvr(pType, ovr.offset);
	    int strengthOvr = strengthMap[polyominoOvr];
	    if (strengthOvr >= minBinding and overlappedOvr == 0) {
		eventSet->addAttachmentEvent(polyominoOvr);
	    }
	}
    }
}
