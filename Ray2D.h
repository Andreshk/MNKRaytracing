#pragma once
#include <cmath>
#include <vector>
#include <algorithm> // std::{min,max}
#include <type_traits> // std::numeric_limits::infinity
#include "..\..\Vassert.h"

bool approx(const float& x, const float& y) { return (std::abs(x - y) < 1e-6f); }

enum class Axis { X, Y };
enum class Space { World, Object };

template <Space S>
struct Point {
    float x, y;
          float& operator[](const Axis a)       { return (a == Axis::X ? x : y); }
    const float& operator[](const Axis a) const { return (a == Axis::X ? x : y); }
};

template <Space S>
struct Vector {
    float x, y;
          float& operator[](const Axis a)       { return (a == Axis::X ? x : y); }
    const float& operator[](const Axis a) const { return (a == Axis::X ? x : y); }
};

template <Space S>
Vector<S> normalize(const Vector<S>& v) {
    const float d = std::sqrt(v.x*v.x + v.y*v.y);
    return { v.x / d, v.y / d };
}

template <Space S>
Point<S> operator+(const Point<S>& p, const Vector<S>& v) { return { p.x + v.x, p.y + v.y }; }
template <Space S>
Point<S>& operator+=(Point<S>& p, const Vector<S>& v) { p.x += v.x; p.y += v.y; return p; }
template <Space S>
Vector<S> operator-(const Point<S>& p0, const Point<S>& p1) { return { p0.x - p1.x, p0.y - p1.y }; }
template <Space S>
Vector<S> operator+(const Vector<S>& v0, const Vector<S>& v1) { return { v0.x + v1.x, v0.y + v1.y }; }
template <Space S>
Vector<S> operator*(const Vector<S>& v, const float k) { return { v.x*k, v.y*k }; }

template <Space S>
struct Ray {
    Point<S> origin;
    Vector<S> dir;
    Ray(const Point<S>& origin, const Vector<S>& dir) : origin{ origin }, dir{ normalize(dir) } {}
};

const float inf = std::numeric_limits<float>::infinity();

struct SegmentHit {
    float t, u;

    SegmentHit() : t{ inf } {} // Infinity denotes an invalid hit (all other t will be smaller)
    SegmentHit(const float t, const float u) : t{ t }, u{ u } {}
    operator bool() const { return (t < inf); }
};

// All segments are in object space
struct Segment {
    Point<Space::Object> p0, p1;

    SegmentHit intersect(const Ray<Space::Object>& ray, const float tmin = 0.f, const float tmax = 1e18f) const {
        const Vector<Space::Object> w = p0 - ray.origin;
        const Vector<Space::Object> segDir = p1 - p0;
        const float denom = (ray.dir.x*segDir.y - ray.dir.y*segDir.x); // same sign as cross-product, i.e. 0 if parallel
        if (approx(denom, 0.f)) { // Ray & segment are parallel
            const bool collinear = approx(w.x*segDir.y - w.y*segDir.x, 0.f);
            if (!collinear) {
                return {}; // Parallel, but different lines -> surely no intersection
            }
            // Edge cases omitted for brevity
            //vassert(false && "to-do");
            return {};
        }
        const float u = (ray.dir.y*w.x - ray.dir.x*w.y) / denom;
        if (u < 0.f || u > 1.f) {
            return {};
        }
        const float t = -(segDir.x*w.y - segDir.y*w.x) / denom;
        return ((t >= tmin && t < tmax) ? SegmentHit{ t,u } : SegmentHit{});
    }
    Vector<Space::Object> geomNormal() const {
        return normalize(Vector<Space::Object>{ p0.y - p1.y, p1.x - p0.x });
    }
};

struct Hit {
    int primIdx;
    float t, u;
    
    Hit() : t{ inf } {} // Infinity denotes an invalid hit (all other t will be smaller)
    Hit(const int primIdx, const float t, const float u) : primIdx{ primIdx }, t{ t }, u{ u } {}
    operator bool() const { return (t < inf); }

    void improveBy(const Hit& other) {
        // If other.t < t, than other.t != inf, but the extra check is for clarity
        if (other && other.t < t) {
            *this = other;
        }
    }
    void improveBy(const int primIdx, const SegmentHit& other) {
        if (other && other.t < t) {
            this->primIdx = primIdx;
            t = other.t;
            u = other.u;
        }
    }
};

template <Space S>
class Box {
    Point<S> pmin, pmax;
public:
    // A default-constructed box is invalid, expansion by a single point makes it valid.
    Box() : pmin{ 1e18f, 1e18f }, pmax{ -1e18f, -1e18f } {}
    // Merge two boxes into one
    Box(const Box& b1, const Box& b2) {
        pmin.x = std::min(b1.pmin.x, b2.pmin.x);
        pmin.y = std::min(b1.pmin.y, b2.pmin.y);
        pmax.x = std::max(b1.pmax.x, b2.pmax.x);
        pmax.y = std::min(b1.pmax.y, b2.pmax.y);
    }
    // Include a given point in a box
    void expand(const Point<S>& p) {
        pmin.x = std::min(pmin.x, p.x);
        pmin.y = std::min(pmin.y, p.y);
        pmax.x = std::max(pmax.x, p.x);
        pmax.y = std::max(pmax.y, p.y);
    }
    // Note: infinity denotes no hit (similar to the Hit structs)
    float intersect(const Ray<S>& ray, float tmin = 0.f, float tmax = 1e18f) const {
        bool hitN = false; // Whether the first wall (relative to ray direction) is hit
        for (const Axis ax : {Axis::X, Axis::Y}) {
            // Note: safe to do divisions by zero; the resulting inf-s will get filtered out by the min/max
            const float t0 = std::min((pmin[ax] - ray.origin[ax]) / ray.dir[ax],
                                      (pmax[ax] - ray.origin[ax]) / ray.dir[ax]);
            const float t1 = std::max((pmin[ax] - ray.origin[ax]) / ray.dir[ax],
                                      (pmax[ax] - ray.origin[ax]) / ray.dir[ax]);
            //tmin = std::max(t0, tmin);
            if (t0 > tmin) {
                tmin = t0;
                hitN = true;
            }
            tmax = std::min(t1, tmax);
            if (tmax <= tmin) {
                return inf;
            }
        }
        return (hitN ? tmin : tmax);
    }
    // Split by a line, perpendicular to the given axis (in 3D we split by a plane, perpendicular to the axis)
    std::pair<Box, Box> split(const Axis ax, const float splitVal) const {
        vassert(pmin[ax] < splitVal && splitVal < pmax[ax]);
        Box left = *this, right = *this;
        left.pmax[ax] = splitVal;
        right.pmin[ax] = splitVal;
        return { left, right };
    }
};

class Mesh {
    std::vector<Point<Space::Object>> vertices; // Set of all vertices
    std::vector<int> indices; // 2 indices per primitive
    std::vector<Vector<Space::Object>> normals; // Set of all normals
    std::vector<int> normalIndices; // 2 indices per primitive
public:
    const Box<Space::Object> box; // Constructed during parsing

    int numPrimitives() const { return indices.size() / 2; }
    Segment getPrimitive(const int idx) const {
        return { vertices[indices[2*idx]],
                 vertices[indices[2*idx+1]] };
    }
    Vector<Space::Object> getNormal(const int idx, const float u) const {
        return normalize(normals[normalIndices[2*idx  ]]*(1-u)
                       + normals[normalIndices[2*idx+1]]*u);
    }
    // For 3D meshes (when index buffers have 3 values per primitive):
    //struct UVCoord { float u, v; }; // Can be Point<Space::UV>
    //std::vector<UVCoord> uvs; // Set of all UVs
    //std::vector<int> ivIndices; // 3 indices per primitive
    //UVCoord Mesh::getUVs(const int idx, const float u, const float v) const {
    //    return uvs[uvIndices[3*idx  ]]*(1-u-v)
    //         + uvs[uvIndices[3*idx+1]]*u
    //         + uvs[uvIndices[3*idx+2]]*v;

    Hit intersect(const Ray<Space::Object>& ray, const float tmin = 0.f, const float tmax = 1e18f) const {
        if (box.intersect(ray, tmin, tmax) != inf) {
            return {};
        }
        Hit bestHit;
        for (int i = 0; i < numPrimitives(); ++i) {
            bestHit.improveBy(i, getPrimitive(i).intersect(ray, tmin, tmax));
            // We can bump tmin and tmax on successful improvement, but there's no perf benefit of doing it
        }
        return bestHit;
    }
};

class Kdtree {
    struct Node { // Note: node size can be reduced drastically
        const bool isLeaf;
        float splitVal;
        int left, right;

        Node(const int left, const int right) // Leaf constructor, from primIndices subrange
            : isLeaf{ true }, left{ left }, right{ right } {}
        Node(const float splitVal) // Non-leaf constructor, from split value (split axis is implicit)
            : isLeaf{ false }, splitVal{ splitVal } {}
    };

    const Mesh& geometry;

    // All tree nodes. Root node is [0], non-leaf nodes contain child node indices in this array.
    // The geometry's bounding box is used to obtain each node's implicit bounding box during traversal.
    std::vector<Node> nodes;
    // Contains reordered primitive indices. Each leaf node points to a range in this array.
    // As an alternative, we can rearrange the geometry's internal index buffer.
    std::vector<int> primIndices;
    // Construction parameters
    static const int maxDepth = 60;
    static const int maxPrimsInLeaf = 15;


    // Invariant: constructs a node, managing the primitives, indexed by the range primIndices[from;to).
    // May rearrange the values in this range. Adds a single node to the array & returns its index.
    int buildRecursive(const int from, const int to, const int depth) {
        if (depth > maxDepth) {
            nodes.push_back({ from, to });
            return (nodes.size() - 1); // The index of the freshly inserted node
        }
        const Axis ax = (depth % 2 ? Axis::Y : Axis::X); // Split axis implicit from depth (not necessary)
        const float splitVal = getMidPoint(from, to, ax);
        const int mid = partition(from, to, ax, splitVal);
        nodes.push_back({ splitVal }); // Insert before the recursive calls (!)
        const int nodeIdx = nodes.size() - 1; // The index of the freshly inserted node
        // Do not make a local Node& node = nodes[nodeIdx] - it'll become dangling after vector reallocations!
        nodes[nodeIdx].left = buildRecursive(from, mid, depth + 1); // Note: always returns nodeIdx+1
        nodes[nodeIdx].right = buildRecursive(mid, to, depth + 1);
        return nodeIdx;
    };

    float getMidPoint(const int from, const int to, const Axis ax) {
        vassert(false && "to-do");
        return -1.f;
    }
    int partition(const int from, const int to, const Axis ax, const float midPoint) {
        // Important: some primitives may end up straddling the split line!
        // For segments it's ok to split them, triangles are better off duplicated.
        // Either way, the indices array may grow, i.e. 'to' should be an int& to account for that!
        vassert(false && "to-do");
        return (from + to) / 2;
    }

    // Invariant: box is the AABB for all primitives in nodes[nodeIdx], and we have already tested whether ray intersects with it.
    Hit intersectRecursive(
        const Ray<Space::Object>& ray, const Box<Space::Object>& box,
        const float tmin, const float tmax, const int nodeIdx, const int depth
    ) const {
        const Node& node = nodes[nodeIdx];
        if (node.isLeaf) {
            Hit bestHit;
            // there would be less confusion if this was an std::range or std::span :)
            for (int i = node.left; i < node.right; ++i) {
                const int primIdx = primIndices[i];
                bestHit.improveBy(primIdx, geometry.getPrimitive(primIdx).intersect(ray, tmin, tmax));
            }
            return bestHit;
        } else {
            const Axis ax = (depth % 2 ? Axis::Y : Axis::X);
            const auto [leftBox, rightBox] = box.split(ax, node.splitVal);
            const bool leftIsNearest = (ray.dir[ax] > 0);
            const int nearIdx = (leftIsNearest ? node.left : node.right);
            const int farIdx = (leftIsNearest ? node.right : node.left);
            const auto& nearBox = (leftIsNearest ? leftBox : rightBox);
            const auto& farBox = (leftIsNearest ? rightBox : leftBox);
            
            // Try intersecting with the nearest box (and on success, its primitives) first.
            // If there is a hit with some primitive inside it, there's no point in intersecting the other box.
            // Important: this is guaranteed by the fact that the two bounding boxes do not overlap!
            // Also, if the ray is entirely before/beyond the split line, we can skip intersecting one of the boxes.
            Hit bestHit;
            if (nearBox.intersect(ray, tmin, tmax) != inf) {
                // Note: it is beneficial to use the box intersection result here to limit the [tmin, tmax] range here.
                // In order to do this, we need to know whether a given box is intersected once or twice (!)
                bestHit.improveBy(intersectRecursive(ray, nearBox, tmin, tmax, nearIdx, depth + 1));
            } else if (farBox.intersect(ray, tmin, tmax) != inf) {
                // Same comment about modifying [tmin;tmax] here, too
                bestHit.improveBy(intersectRecursive(ray, farBox, tmin, tmax, farIdx, depth + 1));
            }
            return bestHit;
        }
    }
public:
    Kdtree(const Mesh& geometry_) : geometry{ geometry_ } {
        const int numPrimitives = geometry.numPrimitives();
        primIndices.resize(numPrimitives);
        // Build our index array (not to be confused with the mesh's one)
        for (int i = 0; i < numPrimitives; ++i) {
            primIndices[i] = i;
        }
        buildRecursive(0, numPrimitives, 0); // This returns 0, i.e. the root node index
    }

    Hit intersect(const Ray<Space::Object>& ray, Hit& hit, const float tmin = 0.f, const float tmax = 1e18f) const {
        // Modifying [tmin;tmax] here might not be needed here
        if (geometry.box.intersect(ray, tmin, tmax) != inf) {
            return intersectRecursive(ray, geometry.box, tmin, tmax, 0, 0);
        } else {
            return {};
        }
    }
};

template <Space S1, Space S2>
struct Transform {
    static_assert(S1 != S2);
    float m[2][2]; // This can be Matrix<S1, S2>
    float off[2];  // This can be Vector<S2>, indicating addition after multiplication
};

// Post-multiplication of matrix & column-vector
template <Space S1, Space S2>
Point<S2> operator*(const Transform<S1, S2>& tm, const Point<S1>& p) {
    // This can be tm.m*p + tm.off;
    return { tm.m[0][0]*p.x + tm.m[0][1]*p.y + tm.off[0],
             tm.m[1][0]*p.x + tm.m[1][1]*p.y + tm.off[1] };
}
template <Space S1, Space S2>
Vector<S2> operator*(const Transform<S1, S2>& tm, const Vector<S1>& p) {
    // This can be just tm.m*p. Note: no translation (!)
    return { tm.m[0][0]*p.x + tm.m[0][1]*p.y,
             tm.m[1][0]*p.x + tm.m[1][1]*p.y };
}

template <Space S1, Space S2>
Transform<S2, S1> inverse(const Transform<S1, S2>& tm) {
    Transform<S2, S1> res;
    const float d = tm.m[0][0]*tm.m[1][1] - tm.m[0][1]*tm.m[1][0];
    res.m[0][0] =  tm.m[1][1]/d;
    res.m[0][1] = -tm.m[0][1]/d;
    res.m[1][0] = -tm.m[1][0]/d;
    res.m[1][1] =  tm.m[0][0]/d;
    // Note: this is the same as multiplying the negative (i.e. inverse) offset by the inverse matrix
    res.off[0] = -(tm.off[0]*res.m[0][0] + tm.off[1]*res.m[0][1]);
    res.off[1] = -(tm.off[0]*res.m[1][0] + tm.off[1]*res.m[1][1]);
    return res;
}

class TransformedMesh {
    const Mesh& geometry; // Note: referenced, not owned

    Box<Space::World> box; // Constructed during parsing, doesn't make Mesh::box obsolete
    Transform<Space::Object, Space::World> tm; // Transforms are also supplied during scene parsing
    Transform<Space::World, Space::Object> itm;
public:
    Hit intersect(const Ray<Space::World>& ray, const float tmin = 0.f, const float tmax = 1e18f) const {
        if (box.intersect(ray, tmin, tmax) != inf) {
            return {};
        }
        const Ray<Space::Object> rayObj{ itm*ray.origin, itm*ray.dir };
        const Hit objHit = geometry.intersect(rayObj, tmin, tmax);
        // Now apply tm to whatever results we extract in object-space
        // (f.e. normals, but carefully). Scalars should not be transformed :)
        return objHit;
    }
};
