#ifndef STOCHASTIC_HPP
#define STOCHASTIC_HPP

#include <climits>
#include <vector>
#include <string>
#include <shared_mutex>
#include <unordered_map>
#include <unordered_set>
#include "polyomino.hpp"
#include "event.hpp"
#include "pugixml.hpp"

struct DimensionRestriction {
    Vec3 min, max;

    DimensionRestriction() {
	min = Vec3(INT_MIN, INT_MIN, INT_MIN);
	max = Vec3(INT_MAX, INT_MAX, INT_MAX);
    }
};

struct ShapeType
{
    std::vector<PolyominoType*> polyominoTypes;
    std::vector<Vec3> coords;

    Vec3 minCoord;
    Vec3 maxCoord;
	
    int id;
};

struct Shape
{
    ShapeType *type;
    Vec3 offset;

    Shape(ShapeType *type=nullptr, Vec3 offset=Vec3(0,0,0))
	: type(type), offset(offset) {}
};

inline bool operator==(const Shape &lhs, const Shape &rhs)
{
    return lhs.type == rhs.type and lhs.offset == rhs.offset;
}

namespace std
{
    template <>
    struct hash<Shape>
    {
	// FNV-1a hash implementation
	std::size_t operator()(const Shape &shape) const noexcept {
	    const unsigned char* tptr = (const unsigned char*)&(shape.type);
	    const unsigned char* xptr = (const unsigned char*)&(shape.offset.x);
	    const unsigned char* yptr = (const unsigned char*)&(shape.offset.x);
	    const unsigned char* zptr = (const unsigned char*)&(shape.offset.x);
	    std::size_t h = FNV_OFFSET;
	    for (unsigned i = 0; i < sizeof(ShapeType*); ++i)
		h = (h ^ tptr[i]) * FNV_PRIME;
	    for (unsigned i = 0; i < sizeof(int); ++i)
		h = (h ^ xptr[i]) * FNV_PRIME;
	    for (unsigned i = 0; i < sizeof(int); ++i)
		h = (h ^ yptr[i]) * FNV_PRIME;
	    for (unsigned i = 0; i < sizeof(int); ++i)
		h = (h ^ zptr[i]) * FNV_PRIME;
	    return h;
	}
    };
}

bool complementaryDomainLabels(const std::string &a, const std::string &b);
Vec3 directionOffset(int direction);
int blockStrength(const Block &b1, const Block &b2, const Vec3 &offset);

struct PrecomputedMaps
{
    std::vector<PolyominoType*> polyominoTypes;
    std::vector<ShapeType*> polyominoShapes;
    std::vector<ShapeType*> shapeTypes;
    std::vector<std::vector<Shape>> neighborLists;
    std::vector<std::vector<Shape>> overlapLists;
    std::vector<std::vector<std::unordered_map<Vec3, int>>> bindingMap;
    std::vector<std::vector<std::pair<PolyominoType*, Vec3>>> potentialBindingList;
    std::unordered_map<Domain, std::vector<std::pair<PolyominoType*, Vec3>>> domainLookup;

    std::vector<std::vector<std::shared_mutex>> bindingMapLocks;

    int getBinding(PolyominoType *pt1, PolyominoType *pt2, const Vec3 &offset);
    
    PrecomputedMaps(unsigned int size) {
	this->bindingMapLocks.resize(size);
	
	for (unsigned int i = 0; i < size; i++) {
	    std::vector<std::shared_mutex> row(size);
	    this->bindingMapLocks[i].swap(row);
	}
    }
};

struct SimParams {
    int threadId;
    int numSteps;
    int maxAttachments;

    double Gmc;
    double Gse;
    double kf;
	
    int bindingThreshold;
    int minBinding;
    DimensionRestriction dimRestrictions;
	
    int outputInterval;
    std::string outputDir;
    std::string outputName;
	
    bool reportNonDeterminism;

    long long rngSeed;
};

class Simulation
{
protected:
    // keeps track of attached polyominoes by shape
    // should only have entries for shapes with attached polyominoes
    std::unordered_map<Shape, Polyomino> locationMap;

    // the following may have entries for unoccupied shapes
    // keeps track of which shapes are overlapped
    std::unordered_map<Shape, int> overlappedMap;
    std::unordered_map<Polyomino, int> strengthMap;
    std::unordered_set<Polyomino> seedSet;

    std::string outputDir;
    std::string outputName;
    
    pugi::xml_document *systemXMLDocument;

    PrecomputedMaps *sharedData;

    int threadId;
    int steps;

    long long rngSeed;

    DimensionRestriction dimRestriction;
    int outputInterval;
    int numSteps;
    

    virtual void initSystemXML() = 0;
    virtual void writeSystemXML(std::string note="") = 0;

public:

    Simulation(
	       PrecomputedMaps *sharedData,
	       const std::unordered_set<Polyomino> &seedPolyominoes,
	       const std::unordered_set<Polyomino> &assemblyPolyominoes,
	       SimParams simParams) :
        dimRestriction(simParams.dimRestrictions),
        outputInterval(simParams.outputInterval), 
        outputDir(simParams.outputDir),
        outputName(simParams.outputName),
        threadId(simParams.threadId),
        rngSeed(simParams.rngSeed),
        numSteps(simParams.numSteps)
    {}
    
    virtual void operator()() = 0;
};

class AbstractSimulation : public Simulation
{
    std::unordered_set<Polyomino> frontier;
    std::vector<Polyomino> frontierList;
    std::unordered_map<Polyomino, int> frontierIds;

    // when frontier locations are added or removed, each grid location will be updated so that non-determinism can be evaluated
    std::unordered_map<Vec3, int> frontierGridLocationCounts;

    int bindingThreshold;
    
    std::mt19937 gen;
    int numSteps;
    bool reportNonDeterminism;

    void addFrontier(Polyomino polyomino);
    void removeFrontier(Polyomino polyomino);

    void checkNonDeterminism();

    void initSystemXML();
    void writeSystemXML(std::string note="");
public:
    AbstractSimulation(PrecomputedMaps *sharedData,
		       const std::unordered_set<Polyomino> &seedPolyominoes,
		       const std::unordered_set<Polyomino> &assemblyPolyominoes,
		       SimParams simParams);

    bool step();
    void addPolyomino(Polyomino polyomino, bool seed=false);

    void operator()();
};

class KineticSimulation : public Simulation
{

    double Gmc, Gse, kf;
    StochasticEventSet *eventSet;

    int minBinding;
    int maxAttachments;

    double time;

    void seedRng(int seed);

    void initSystemXML();
    void writeSystemXML(std::string note="");
    
public:

    KineticSimulation(PrecomputedMaps *sharedData,
		      const std::unordered_set<Polyomino> &seedPolyominoes,
		      const std::unordered_set<Polyomino> &assemblyPolyominoes,
		      SimParams simParams);

    void step();
    void addPolyomino(Polyomino polyomino, bool seed=false);
    void removePolyomino(Polyomino polyomino);

    void operator()();
};

#endif
