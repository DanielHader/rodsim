#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string.h>
#include <unordered_set>
#include <utility>
#include "assembly.hpp"
#include <thread>

Assembly::Assembly(const std::string &filename)
{
    this->simulationType = SIM_KINETIC;
    this->reportNonDeterminism = false;
    this->rngSeed = -1;

	this->minBinding = -1;
	int temp_min_binding = -1;
    this->parseConfFile(filename);
	if (this->minBinding > -1) {
		temp_min_binding = this->minBinding;
	}
	
	// use the min_binding value in the config file over the value in the system file
    this->parseSystemFile(this->confPath + this->system);
	if (temp_min_binding != -1) {
		this->minBinding = temp_min_binding;
	}

    std::cout << "found " << this->polyominoTypes.size() << " polyomino types" << std::endl << std::endl;

    std::cout << "baking precomputation maps from type data" << std::endl;
    this->computeShapes();
    this->computeBoundingBoxes();
    this->computeNeighborLists();
    this->computeDomainLookup();
    this->computePotentialBindingLists();
    this->initializeBindingMap();
    std::cout << std::endl;

    this->sharedData = new PrecomputedMaps(this->polyominoTypes.size());

    // move precomputed maps to shared data
    this->sharedData->polyominoTypes = std::move(this->polyominoTypes);
    this->sharedData->polyominoShapes = std::move(this->polyominoShapes);
    this->sharedData->shapeTypes = std::move(this->shapeTypes);
    this->sharedData->neighborLists = std::move(this->neighborLists);
    this->sharedData->overlapLists = std::move(this->overlapLists);
    this->sharedData->bindingMap = std::move(this->bindingMap);
    this->sharedData->potentialBindingList = std::move(this->potentialBindingList);
    this->sharedData->domainLookup = std::move(this->domainLookup);
}

void Assembly::run()
{
    // initialize ensemble
    std::cout << "initializing ensemble threads and beginning simulations" << std::endl;
    // std::vector<boost::thread*> threads;

    std::vector<std::thread> threads;
    std::string filename = confPath;
    if (outputDir != "")
	filename += outputDir + "/";

    if (this->simulationType == SIM_KINETIC)
    {
		SimParams simParams;
		simParams.Gmc = this->Gmc;
		simParams.Gse = this->Gse;
		simParams.kf = this->kf;
		simParams.dimRestrictions = this->dimRestriction;
		simParams.minBinding = this->minBinding;
		simParams.outputInterval = this->outputInterval;
		simParams.outputDir = filename;
		simParams.outputName = this->outputName;
		simParams.numSteps = this->numSteps;
		simParams.maxAttachments = this->maxAttachments;
		simParams.rngSeed = this->rngSeed;

		std::vector<KineticSimulation> sims;
		for (int i = 0; i < this->ensembleSize; i++)
		{
			simParams.threadId = i;

			sims.push_back(KineticSimulation(this->sharedData, this->seedPolyominoes, this->assemblyPolyominoes, simParams));
			threads.push_back(std::thread(sims[i]));
		}
    } else {
		SimParams simParams;
		simParams.bindingThreshold = this->bindingThreshold;
		simParams.dimRestrictions = this->dimRestriction;
		simParams.minBinding = this->minBinding;
		simParams.outputInterval = this->outputInterval;
		simParams.outputDir = filename;
		simParams.outputName = this->outputName;
		simParams.numSteps = this->numSteps;
		simParams.reportNonDeterminism = this->reportNonDeterminism;
		simParams.rngSeed = this->rngSeed;

		std::vector<AbstractSimulation> sims;
		for (int i = 0; i < this->ensembleSize; i++)
		{
			simParams.threadId = i;

			sims.push_back(AbstractSimulation(this->sharedData, this->seedPolyominoes, this->assemblyPolyominoes, simParams));

			threads.push_back(std::thread(sims[i]));
		}
    }

    for (std::thread &thread : threads)
	thread.join();
}

void Assembly::parseConfFile(const std::string &filename)
{
    std::ifstream confStream(filename);
	
    if (!confStream.is_open())
	throw std::runtime_error("unable to open conf file: " + filename);

    std::cout << "parsing conf file: " << filename << std::endl;
	
    this->confPath = "";
    std::size_t pos = filename.find_last_of("/");
    if (pos != std::string::npos)
		this->confPath = filename.substr(0, pos+1);
	
    std::string line;
    std::string param;
    std::string value;
    while (std::getline(confStream, line)) {
		pos = line.find("=");

		if (pos != std::string::npos) {
			param = line.substr(0, pos);
			value = line.substr(pos + 1);
			this->parseConfParam(param, value);
		}
    }

    confStream.close();
}

int Assembly::parseConfParam(const std::string &param, const std::string &value) {
    if (param == "num_steps") {
		this->numSteps = std::stoi(value);
    } else if (param == "max_attachments") {
		this->maxAttachments = std::stoi(value);
    } else if (param == "min_binding") {
	this->minBinding = std::stoi(value);
    } else if (param == "output_interval") {
		this->outputInterval = std::stoi(value);
    } else if (param == "system") {
		this->system = value;
    } else if (param == "output_name") {
		this->outputName = value;
    } else if (param == "output_directory") {
		this->outputDir = value;
    } else if (param == "ensemble_size") {
		this->ensembleSize = std::stoi(value);
    } else if (param == "sim_type") {
	if (value == "abstract")
	    this->simulationType = SIM_ABSTRACT;
	else
	    this->simulationType = SIM_KINETIC;
    } else if (param == "report_nondeterminism") {
		if (value == "true" or value == "True") {
			this->reportNonDeterminism = true;
		}
    } else if (param == "rng_seed") {
		this->rngSeed = std::stoll(value);
    } else if (param == "runtime_memoization") {
		if (value == "true" or value == "True") {
			this->runtimeMemoization = true;
		}
    } else {
		return 1;
    }
	
    return 0;
}

void Assembly::writeSystemFile(const std::string &filename) {
    pugi::xml_node sysNode = systemXMLDocument.child("PolyominoSystem");
    /*
    sysNode.remove_child("time");
    pugi::xml_node timeNode = sysNode.prepend_child("time");
    timeNode.append_child(pugi::node_pcdata).set_value(std::to_string(sim->time).c_str());

    sysNode.remove_child("steps");
    pugi::xml_node stepsNode = sysNode.prepend_child("steps");
    stepsNode.append_child(pugi::node_pcdata).set_value(std::to_string(sim->steps).c_str());
    */
    sysNode.remove_child(systemXMLDocument.child("assembly"));
    pugi::xml_node assemblyNode = sysNode.append_child("assembly");
    pugi::xml_node polyominoesNode = assemblyNode.append_child("Polyominoes");

    /*
    for (auto kv : sim->locationMap) {
	Polyomino p = kv.second;

	pugi::xml_node polyominoNode = polyominoesNode.append_child("Polyomino");
	pugi::xml_node typeNode = polyominoNode.append_child("PolyominoType");
	typeNode.append_child(pugi::node_pcdata).set_value(p.type->name.c_str());
	pugi::xml_node translationNode = polyominoNode.append_child("translation");
	translationNode.append_child(pugi::node_pcdata).set_value(p.offset.toString().c_str());
    }
    */

    systemXMLDocument.save_file(filename.c_str());
}

void Assembly::parseSystemFile(const std::string &filename)
{
    pugi::xml_parse_result result;
    result = this->systemXMLDocument.load_file(filename.c_str());

    if (!result)
	throw std::runtime_error("unable to load XML file: " + filename);

    std::cout << "parsing system file: " << filename << std::endl;
	
    pugi::xml_node sysNode = systemXMLDocument.child("PolyominoSystem");

    this->Gmc = sysNode.child("Gmc").text().as_double(19.1);
    this->Gse = sysNode.child("Gse").text().as_double(10.0);

    if (sysNode.child("forwardRate").empty())
		this->kf = sysNode.child("forward_rate").text().as_double(7.88e5);
    else
		this->kf = sysNode.child("forwardRate").text().as_double(7.88e5);

    this->bindingThreshold = sysNode.child("bindingThreshold").text().as_int(2);
    this->minBinding = sysNode.child("min_binding").text().as_int(0);
    this->rotation = sysNode.child("rotation").text().as_bool(false);

    pugi::xml_node typesNode = sysNode.child("PolyominoTypes");
    for (pugi::xml_node typeNode : typesNode.children("PolyominoType")) {
		this->parsePolyominoType(typeNode);
    }

	pugi::xml_node dimRestrictionsNode = sysNode.child("dim_restrictions");
	std::string dimRestrictionsString = dimRestrictionsNode.text().get();
	this->parseDimRestrictions(dimRestrictionsString);
	
    pugi::xml_node seedNode = sysNode.child("seed");
    pugi::xml_node psNode = seedNode.child("Polyominoes");
    for (pugi::xml_node polyNode : psNode.children("Polyomino")) {
		
	std::string name = polyNode.child("PolyominoType").text().get();
	if (this->polyominoTypesByName.count(name) == 0)
	    throw std::runtime_error("no polyomino type with name: " + name);
	PolyominoType *polyominoType = this->polyominoTypesByName[name];

	const char *coord_str = polyNode.child("translation").text().get();
	Vec3 offset = this->parseCoord(coord_str);

	Polyomino polyomino(polyominoType, offset);
		this->seedPolyominoes.insert(polyomino);
    }

    pugi::xml_node assemblyNode = sysNode.child("assembly");
    pugi::xml_node paNode = assemblyNode.child("Polyominoes");
    for (pugi::xml_node polyNode : paNode.children("Polyomino")) {

		std::string name = polyNode.child("PolyominoType").text().get();
		if (this->polyominoTypesByName.count(name) == 0)
			throw std::runtime_error("no polyomino type with name: " + name);
		PolyominoType *polyominoType = this->polyominoTypesByName[name];
			
		const char *coord_str = polyNode.child("translation").text().get();
		Vec3 offset = this->parseCoord(coord_str);

		Polyomino polyomino(polyominoType, offset);
		if (!seedPolyominoes.count(polyomino))
		{
			this->assemblyPolyominoes.insert(polyomino);
		}
    }
}

void Assembly::parseDimRestrictions(const std::string &s) {
	const std::regex dimRegex("\\(\\s*(None|\\[(\\-?\\d+),\\s*(\\-?\\d+)\\])\\s*,\\s*(None|\\[(\\-?\\d+),\\s*(\\-?\\d+)\\])\\s*,\\s*(None|\\[(\\-?\\d+),\\s*(\\-?\\d+)\\])\\s*\\)");
	
	std::smatch match;

	// TODO: perform validation on numbers
	if (std::regex_match(s, match, dimRegex)) {
		if (match[1].str() != "None") {
			try {
				this->dimRestriction.min.x = std::stoi(match[2].str());
				this->dimRestriction.max.x = std::stoi(match[3].str());
			} catch (std::exception err) {
				std::cout << "failed to read dim_restrictions x: " << match[1].str() << std::endl;
			}
		}
		if (match[4].str() != "None") {
			try {
				this->dimRestriction.min.y = std::stoi(match[5].str());
				this->dimRestriction.max.y = std::stoi(match[6].str());
			} catch (std::exception err) {
				std::cout << "failed to read dim_restrictions y: " << match[4].str() << std::endl;
			}
		}
		if (match[7].str() != "None") {
			try {
				this->dimRestriction.min.z = std::stoi(match[8].str());
				this->dimRestriction.max.z = std::stoi(match[9].str());
			} catch (std::exception err) {
				std::cout << "failed to read dim_restrictions z: " << match[4].str() << std::endl;
			}
		}

		std::cout << this->dimRestriction.min.toString() << std::endl;
		std::cout << this->dimRestriction.max.toString() << std::endl;
	}
}

void Assembly::parsePolyominoType(pugi::xml_node typeNode)
{
    PolyominoType *type = new PolyominoType;
    type->name = typeNode.child("name").text().as_string();
    type->color = typeNode.child("color").text().as_string();
    type->concentration = typeNode.child("concentration").text().as_double(1.0);

    pugi::xml_node blocksNode = typeNode.child("blocks");
    for (pugi::xml_node blockNode : blocksNode.children("block")) {
		Block block;
			
		const char *coord_str = blockNode.child("coords").text().get();
		block.coord = this->parseCoord(coord_str);

		pugi::xml_node domainsNode = blockNode.child("domains");
		for (pugi::xml_node domainNode : domainsNode.children("domain")) {
			Domain domain;
				
			domain.label = domainNode.child("label").text().as_string();
			domain.strength = domainNode.child("strength").text().as_int(1);

			const char *dir_str = domainNode.child("direction").text().get();
			domain.direction = this->parseDirection(dir_str);

			block.domains.push_back(domain);
		}
		type->blocks.push_back(block);
    }

    type->id = this->polyominoTypes.size();
    this->polyominoTypes.push_back(type);
    if (this->polyominoTypesByName.count(type->name))
		throw std::runtime_error("duplicate polyomino type name: " + type->name);
	
    this->polyominoTypesByName[type->name] = type;
}

Vec3 Assembly::parseCoord(const char *coord_str)
{
    Vec3 coord;
    if (sscanf(coord_str, "(%d,%d,%d)", &coord.x, &coord.y, &coord.z) == 3) {
		return coord;
    } else if (sscanf(coord_str, "(%d,%d)", &coord.x, &coord.y) == 2) {
		coord.z = 0;
		return coord;
    } else {
		throw std::runtime_error(std::string("invalid coord string: ")+coord_str);
    }
}

int Assembly::parseDirection(const char *dir_str)
{
    if (strlen(dir_str) < 1)
	throw std::runtime_error(std::string("invalid direction string: ")+dir_str);
	
    switch (dir_str[0]) {
    case 'N': case 'n': return NORTH;
    case 'E': case 'e': return EAST;
    case 'S': case 's': return SOUTH;
    case 'W': case 'w': return WEST;
    case 'U': case 'u': return UP;
    case 'D': case 'd': return DOWN;
    default:
	throw std::runtime_error(std::string("invalid direction string: ")+dir_str);
    }
}

// precomputation functions
bool isMatchingShape(ShapeType *sType, PolyominoType *pType)
{
    if (sType->coords.size() != pType->blocks.size())
	return false;
	
    for (const Block &block : pType->blocks) {
	bool found = false;
	for (const Vec3 &coord : sType->coords) {
	    if (block.coord == coord) {
		found = true;
		break;
	    }
	}
	if (!found)
	    return false;
    }
    return true;
}

inline void coordMin(Vec3 &a, const Vec3 &b) {
    a.x = std::min(a.x, b.x);
    a.y = std::min(a.y, b.y);
    a.z = std::min(a.z, b.z);
}

inline void coordMax(Vec3 &a, const Vec3 &b) {
    a.x = std::max(a.x, b.x);
    a.y = std::max(a.y, b.y);
    a.z = std::max(a.z, b.z);
}

void Assembly::computeShapes()
{
    std::cout << "  computing unique polyomino shapes... " << std::flush;

    std::vector<Vec3> polyominoTranslations;

    // translate all shapes so that they're lowest coordinate is (0,0,0)
    for (PolyominoType *pType : this->polyominoTypes) {
	Vec3 minCoord = pType->blocks[0].coord;

	for (const Block &block : pType->blocks) {
	    coordMin(minCoord, block.coord);
	}

	polyominoTranslations.push_back(minCoord);

	for (Block &block : pType->blocks) {
	    block.coord -= minCoord;
	}
    }

    std::unordered_set<Polyomino> seed = std::unordered_set<Polyomino>(std::move(this->seedPolyominoes));
    std::unordered_set<Polyomino> assembly = std::unordered_set<Polyomino>(std::move(this->assemblyPolyominoes));

    this->seedPolyominoes.clear();
    this->assemblyPolyominoes.clear();

    for (Polyomino p : seed) {
	p.offset += polyominoTranslations[p.type->id];
	this->seedPolyominoes.insert(p);
    }

    for (Polyomino p : assembly) {
	p.offset += polyominoTranslations[p.type->id];
	this->assemblyPolyominoes.insert(p);
    }

    for (PolyominoType *pType : this->polyominoTypes) {
	ShapeType *foundType = nullptr;
	for (ShapeType *sType : this->shapeTypes) {
	    if (isMatchingShape(sType, pType)) {
		foundType = sType;
		sType->polyominoTypes.push_back(pType);
	    }
	}

	if (foundType == nullptr) {
	    foundType = new ShapeType();
	    foundType->polyominoTypes.push_back(pType);
	    foundType->id = this->shapeTypes.size();
	    for (const Block &block : pType->blocks)
		foundType->coords.push_back(block.coord);

	    this->shapeTypes.push_back(foundType);
	}

	polyominoShapes.push_back(foundType);
    }



    std::cout << "done (" << shapeTypes.size() << " unique shapes)" << std::endl;
}

void Assembly::computeBoundingBoxes() {
    std::cout << "  computing shape bounding boxes... " << std::flush;

    for (ShapeType *sType : this->shapeTypes) {
	sType->minCoord = sType->coords.at(0);
	sType->maxCoord = sType->coords.at(0);
		
	for (const Vec3 &coord : sType->coords) {
	    coordMin(sType->minCoord, coord);
	    coordMax(sType->maxCoord, coord);
	}


    }

    std::cout << "done" << std::endl;
}

inline int taxicab_dist(const Vec3 &a, const Vec3 &b) {
    return std::abs(a.x - b.x) + std::abs(a.y - b.y) + std::abs(a.z - b.z);
}

// returns 2 if overlaps, 1 if neighbors, 0 if neither
int determine_overlap(const ShapeType *s1, const Shape &s2) {
    int ret = 0;
    for (const Vec3 &coord2 : s2.type->coords) {
	Vec3 coord2_off = coord2 + s2.offset;
	for (const Vec3 &coord1 : s1->coords) {
	    int dist = taxicab_dist(coord1, coord2_off);
	    if (dist == 0)
		return 2;
	    else if (dist == 1)
		ret = 1;
	}
    }
    return ret;
}

void Assembly::checkNeighbors(ShapeType *st1, ShapeType *st2)
{
    Vec3 minDiff = st1->minCoord - st2->maxCoord - Vec3(1, 1, 1);
    Vec3 maxDiff = st1->maxCoord - st2->minCoord + Vec3(1, 1, 1);
	
    for (int x = minDiff.x; x <= maxDiff.x; ++x) {
	for (int y = minDiff.y; y <= maxDiff.y; ++y) {
	    for (int z = minDiff.z; z <= maxDiff.z; ++z) {
		Shape s(st2, Vec3(x, y, z));
		int cmp = determine_overlap(st1, s);
		if (cmp == 2) {
		    this->overlapLists.at(st1->id).push_back(s);
		} else if (cmp == 1) {
		    this->neighborLists.at(st1->id).push_back(s);
		}
	    }
	}
    }
}

void Assembly::computeNeighborLists()
{
    std::cout << "  computing shape overlap and neighbor offsets... " << std::flush;
    this->overlapLists.resize(this->shapeTypes.size());
    this->neighborLists.resize(this->shapeTypes.size());
    for (ShapeType *st1 : this->shapeTypes) {
	for (ShapeType *st2 : this->shapeTypes) {
	    this->checkNeighbors(st1, st2);
	}
    }
    std::cout << "done" << std::endl;
}

int relativeStrength(PolyominoType *pt1, PolyominoType *pt2, const Vec3 &offset)
{
    int strength = 0;
    for (const Block &b1 : pt1->blocks) {
	for (const Block &b2 : pt2->blocks) {
	    strength += blockStrength(b1, b2, offset);
	}
    }

    return strength;
}

void Assembly::computeDomainLookup() {
    std::cout << "  computing lookup table for polyomino types by domain... " << std::flush;

    for (PolyominoType *pt : this->polyominoTypes) {
	for (Block block : pt->blocks) {
	    Vec3 coord = block.coord;
	    for (Domain domain : block.domains) {
		if (this->domainLookup.count(domain) == 0) {
		    this->domainLookup[domain] = std::vector<std::pair<PolyominoType*, Vec3>>();
		}

		this->domainLookup[domain].push_back(std::make_pair(pt, coord));
	    }
	}
    }

    std::cout << "done" << std::endl;
}

Domain getComplementaryDomain(Domain domain) {
    Domain comp;
    
    if (domain.label[domain.label.length() - 1] != '*') {
	comp.label = domain.label + "*";
    } else {
	comp.label = domain.label.substr(0, domain.label.length() - 1);
    }


    comp.strength = domain.strength;
    comp.direction = OPP_DIR(domain.direction);

    return comp;
}

void Assembly::computePotentialBindingLists() {
    std::cout << "  computing potential binding pairs... " << std::flush;

    this->potentialBindingList.resize(this->polyominoTypes.size());

    /*
    int strength;
    for (PolyominoType *pt1 : this->polyominoTypes) {
	ShapeType *st1 = this->polyominoShapes[pt1->id];
	for (Shape s2 : this->neighborLists.at(st1->id)) {
	    Vec3 offset = s2.offset;
	    ShapeType *st2 = s2.type;
	    for (PolyominoType *pt2 : st2->polyominoTypes) {
		strength = relativeStrength(pt1, pt2, offset);
		if (strength > 0) {
		    this->potentialBindingList[pt1->id].push_back(std::make_pair(pt2, offset));
		}
	    }
	}
    }
    */

    int strength;
    for (PolyominoType *pt1 : this->polyominoTypes) {

	for (Block block : pt1->blocks) {
	    Vec3 coord = block.coord;
	    for (Domain domain : block.domains) {
		Vec3 dirOffset = directionOffset(domain.direction);
		
		Domain comp = getComplementaryDomain(domain);
		if (this->domainLookup.count(comp) == 0) continue;
		
		for (auto pair : this->domainLookup[comp]) {
		    PolyominoType *pt2 = pair.first;
		    Vec3 offset = pair.second;
		    strength = domain.strength;

		    this->potentialBindingList[pt1->id].push_back(std::make_pair(pt2, coord + dirOffset - offset));
		    //strength = relativeStrength(pt1, pt2, nbrCoord - offset);
		    
		}
	    }
	}
    }

    std::cout << "done" << std::endl;
}

void Assembly::computeBindingMap() {
    
    std::cout << "  computing binding strength of neighboring shapes... " << std::flush;

    this->bindingMap.resize(this->polyominoTypes.size());
    for (unsigned i = 0; i < this->polyominoTypes.size(); ++i)
	this->bindingMap.at(i).resize(this->polyominoTypes.size());

    int strength;
    for (PolyominoType *pt1 : this->polyominoTypes) {
	ShapeType *st1 = this->polyominoShapes[pt1->id];
	for (Shape s2 : this->neighborLists.at(st1->id)) {
	    Vec3 offset = s2.offset;
	    ShapeType *st2 = s2.type;
	    for (PolyominoType *pt2 : st2->polyominoTypes) {
		strength = relativeStrength(pt1, pt2, offset);
		bindingMap.at(pt1->id).at(pt2->id)[offset] = strength;
	    }
	}
    }
    
    std::cout << "done" << std::endl;
}

void Assembly::initializeBindingMap() {

    std::cout << "  initializing binding memo... " << std::flush;

    this->bindingMap.resize(this->polyominoTypes.size());
    for (unsigned i = 0; i < this->polyominoTypes.size(); ++i) {
	this->bindingMap.at(i).resize(this->polyominoTypes.size());
    }
    
    std::cout << "done" << std::endl;
}

int Assembly::memoizedComputeBinding(PolyominoType *pt1, PolyominoType *pt2, Vec3 offset) {
    if (bindingMap.at(pt1->id).at(pt2->id).count(offset) == 0) {
        int strength = relativeStrength(pt1, pt2, offset);
        bindingMap.at(pt1->id).at(pt2->id)[offset] = strength;
    }
    return bindingMap.at(pt1->id).at(pt2->id)[offset];
}
