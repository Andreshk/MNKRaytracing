#pragma once
#include <cmath>
#include <type_traits>
#include "..\..\Vassert.h"

namespace Ray3D {
bool approx(const float& x, const float& y) { return (std::abs(x - y) < 1e-6f); }

struct Point {
    float x, y, z;
};

struct Vector {
    float x, y, z;
};

Point operator+(const Point& p, const Vector& v) { return { p.x + v.x, p.y + v.y, p.z + v.z }; }
Point operator-(const Point& p, const Vector& v) { return { p.x - v.x, p.y - v.y, p.z - v.z }; }
Vector operator-(const Point& p0, const Point& p1) { return { p0.x - p1.x, p0.y - p1.y, p0.z - p1.z }; }
Vector operator*(const float k, const Vector& v) { return { k*v.x, k*v.y, k*v.z }; }
float dot(const Vector& v0, const Vector& v1) {
    return v0.x*v1.x + v0.y*v1.y + v0.z*v1.z;
}
Vector normalize(const Vector& v) {
    const float d = std::sqrt(dot(v, v)); // Also vector length
    return { v.x / d, v.y / d, v.z / d };
}
Vector cross(const Vector& v0, const Vector& v1) {
    return { (v0.y*v1.z - v0.x*v1.y), (v0.z*v1.x - v0.x*v1.z), (v0.x*v1.y - v0.y*v1.x) };
}

struct Ray {
    Point origin;
    Vector dir;
    Ray(const Point& origin, const Vector& dir) : origin{ origin }, dir{ normalize(dir) } {}
};

const float inf = std::numeric_limits<float>::infinity();

struct TriangleHit {
    float t, u, v;

    TriangleHit() : t{ inf } {} // Infinity denotes an invalid hit (all other t will be smaller)
    TriangleHit(const float t, const float u, const float v) : t{ t }, u{ u }, v{ v } {}
    operator bool() const { return (t < inf); }
};

struct Triangle {
    Point v0;
    Vector e1, e2; // Note: can be vertices

    TriangleHit intersect(const Ray& ray, const float tmin = 0.f, const float tmax = 1e18f) const {
        // First intersect the ray with the triangle's plane
        const Vector n = normalize(cross(e1, e2));
        const float denom = dot(n, ray.dir);
        if (approx(denom, 0.f)) { // Ray & plane are collinear or parallel
            // Boring edge cases here, assume no intersection
            vassert(false);
            return {};
        }
        const float t = dot(n, v0 - ray.origin) / denom;
        if (t < tmin || t >= tmax) {
            return {};
        }
        // The intersection point
        const Point inters = ray.origin + t*ray.dir;
        const float e1e2 = dot(e1, e2);
        const float e1e1 = dot(e1, e1);
        const float e2e2 = dot(e2, e2);
        const float denom2 = e1e2*e1e2 - e1e1*e2e2;
        vassert(!approx(denom2, 0.f)); // Triangle is degenerate
        const Vector w = inters - v0;
        const float u = (e1e2*dot(w, e2) - e2e2*dot(w, e1)) / denom2;
        const float v = (e1e2*dot(w, e1) - e1e1*dot(w, e2)) / denom2;
        return ((u >= 0 && v >= 0 && u + v <= 1) ? TriangleHit{ t, u, v } : TriangleHit{});
    }
    Vector geomNormal() const {
        return normalize(cross(e1, e2));
    }
};
} // namespace Ray3D
