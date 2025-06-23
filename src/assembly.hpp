#ifndef ASSEMBLY_HPP
#define ASSEMBLY_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "polyomino.hpp"
#include "stochastic.hpp"
#include "pugixml.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#define SIM_ABSTRACT 1
#define SIM_KINETIC 0

struct Assembly {
    //StochasticSim *sim;
    int simulationType;
	
    // all polyomino types
    std::vector<PolyominoType*> polyominoTypes;
    std::unordered_map<std::string, PolyominoType*> polyominoTypesByName;
	
    std::unordered_set<Polyomino> seedPolyominoes;
    std::unordered_set<Polyomino> assemblyPolyominoes;

    double Gmc;
    double Gse;
    double kf;

    DimensionRestriction dimRestriction;
    int minBinding;
    int bindingThreshold;

    bool rotation;
    bool runtimeMemoization;
	
    long long rngSeed;

    int numSteps;
    int maxAttachments;
    int outputInterval;
    int ensembleSize;

    std::string confPath;
    std::string system;
    std::string outputName;
    std::string outputDir;

    bool reportNonDeterminism;

    PrecomputedMaps *sharedData;

    // precomputation maps
    std::vector<ShapeType*> polyominoShapes;
    std::vector<ShapeType*> shapeTypes;
    std::vector<std::vector<Shape>> neighborLists;
    std::vector<std::vector<Shape>> overlapLists;
    
    std::vector<std::vector<std::unordered_map<Vec3, int>>> bindingMap;
    std::vector<std::vector<std::pair<PolyominoType*, Vec3>>> potentialBindingList;
    std::unordered_map<Domain, std::vector<std::pair<PolyominoType*, Vec3>>> domainLookup;
    
    pugi::xml_document systemXMLDocument;
	
    Assembly(const std::string &filename);
    ~Assembly();

    void run();
private:
    int parseConfParam(const std::string &param, const std::string &value);
    void parseConfFile(const std::string &filename);
    void parseSystemFile(const std::string &filename);
    void writeSystemFile(const std::string &filename);

    void parsePolyominoType(pugi::xml_node);
    Vec3 parseCoord(const char*);
    int parseDirection(const char*);

    // functions for precomputation
    void computeShapes();
    void computeBoundingBoxes();
    void checkNeighbors(ShapeType *st1, ShapeType *st2);
    void computeNeighborLists();
    void computePotentialBindingLists();
    void computeBindingMap();
    void computeDomainLookup();

    void initializeBindingMap();
    int memoizedComputeBinding(PolyominoType *pt1, PolyominoType *pt2, Vec3 offset);

};

#endif
