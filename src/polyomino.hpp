#ifndef POLYOMINO_HPP
#define POLYOMINO_HPP

#include <vector>
#include <string>
#include <utility>
#include <functional>
// #include <boost/functional/hash.hpp>

#define NORTH 0
#define SOUTH 5
#define EAST  1
#define WEST  4
#define UP    2
#define DOWN  3

#define OPP_DIR(d) (5-d)

struct Vec3
{
    int x, y, z;
    Vec3(int x=0, int y=0, int z=0) : x(x), y(y), z(z) {}

    inline std::string toString() const {
	return "(" + std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z) + ")";
    }
};

inline Vec3 directionOffset(int direction)
{
    switch (direction)
    {
    case NORTH : return Vec3( 0, 1, 0);
    case SOUTH : return Vec3( 0,-1, 0);
    case EAST  : return Vec3( 1, 0, 0);
    case WEST  : return Vec3(-1, 0, 0);
    case UP    : return Vec3( 0, 0, 1);
    case DOWN  : return Vec3( 0, 0,-1);
    default    : return Vec3( 0, 0, 0);
    }
}

struct Domain
{
    std::string label;
    int direction;
    int strength;
};

struct Block
{
    Vec3 coord;
    std::vector<Domain> domains;
};

struct PolyominoType
{
    std::string name;
    std::string color;

    double concentration;

    std::vector<Block> blocks;

    int id;
};

struct Polyomino
{
    PolyominoType *type;
    Vec3 offset;

    int initialBindStrength;
    int currentBindStrength;

    int stepAdded;
    double timeAdded;

    Polyomino(PolyominoType *type=nullptr, Vec3 offset=Vec3(0,0,0))
	: type(type), offset(offset) {}
};

inline bool operator==(const Vec3 &lhs, const Vec3 &rhs)
{ return lhs.x == rhs.x and lhs.y == rhs.y and lhs.z == rhs.z; }

inline Vec3 operator+(const Vec3 &lhs, const Vec3 &rhs)
{ return Vec3(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z); }

inline Vec3 operator-(const Vec3 &lhs, const Vec3 &rhs)
{ return Vec3(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z); }

inline Vec3& operator+=(Vec3 &lhs, const Vec3 &rhs)
{
    lhs.x += rhs.x;
    lhs.y += rhs.y;
    lhs.z += rhs.z;
    return lhs;
}

inline Vec3& operator-=(Vec3 &lhs, const Vec3 &rhs)
{
    lhs.x -= rhs.x;
    lhs.y -= rhs.y;
    lhs.z -= rhs.z;
    return lhs;
}

inline bool operator>=(const Vec3 &lhs, const Vec3 &rhs)
{
    return lhs.x >= rhs.x and lhs.y >= rhs.y and lhs.z >= rhs.z;
}

inline bool operator<=(const Vec3 &lhs, const Vec3 &rhs)
{
    return lhs.x <= rhs.x and lhs.y <= rhs.y and lhs.z <= rhs.z;
}

inline bool operator>(const Vec3 &lhs, const Vec3 &rhs)
{
    return lhs.x > rhs.x and lhs.y > rhs.y and lhs.z > rhs.z;
}

inline bool operator<(const Vec3 &lhs, const Vec3 &rhs)
{
    return lhs.x < rhs.x and lhs.y < rhs.y and lhs.z < rhs.z;
}

inline bool operator==(const Domain &lhs, const Domain &rhs) {
    return lhs.label == rhs.label and lhs.direction == rhs.direction and lhs.strength == rhs.strength;
}

inline bool operator==(const Polyomino &lhs, const Polyomino &rhs)
{ return lhs.type == rhs.type and lhs.offset == rhs.offset; }

constexpr std::size_t FNV_PRIME = sizeof(std::size_t) == 8 ? 1099511628211u : 16777619u;
constexpr std::size_t FNV_OFFSET = sizeof(std::size_t) == 8 ? 14695981039346656037u : 2166136261u;

inline void hash_combine(std::size_t &s, const std::size_t &h)
{
    s ^= h + 0x9e3779b9 + (s << 6) + (s >> 2);
}

namespace std
{
    template <>
    struct hash<Vec3>
    {
	// FNV-1a hash
	std::size_t operator()(const Vec3 &vec) const noexcept {
	    std::size_t x_hash = std::hash<int>{}(vec.x);
	    std::size_t y_hash = std::hash<int>{}(vec.y);
	    std::size_t z_hash = std::hash<int>{}(vec.z);

	    std::size_t seed = 0;
	    hash_combine(seed, x_hash);
	    hash_combine(seed, y_hash);
	    hash_combine(seed, z_hash);

	    return seed;
	}
    };

    template <>
    struct hash<Domain>
    {
	// FNV-1a hash
	std::size_t operator()(const Domain &domain) const noexcept {
	    std::size_t label_hash = std::hash<std::string>{}(domain.label);
	    std::size_t direction_hash = std::hash<int>{}(domain.direction);
	    std::size_t strength_hash = std::hash<int>{}(domain.strength);

	    std::size_t seed = 0;
	    hash_combine(seed, label_hash);
	    hash_combine(seed, direction_hash);
	    hash_combine(seed, strength_hash);

	    return seed;
	}
    };

    template <>
    struct hash<Polyomino>
    {
	// FNV-1a hash
	std::size_t operator()(const Polyomino &polyomino) const noexcept {
	    std::size_t type_hash = std::hash<int>{}(polyomino.type->id);
	    std::size_t x_hash = std::hash<int>{}(polyomino.offset.x);
	    std::size_t y_hash = std::hash<int>{}(polyomino.offset.y);
	    std::size_t z_hash = std::hash<int>{}(polyomino.offset.z);
			
	    std::size_t seed = 0;
	    hash_combine(seed, type_hash);
	    hash_combine(seed, x_hash);
	    hash_combine(seed, y_hash);
	    hash_combine(seed, z_hash);

	    return seed;
	}
    };
}

#endif
 
