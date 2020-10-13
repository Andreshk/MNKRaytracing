#include <iostream>
#include "Ray2D.h"
#include "Ray3D.h"

// Incomplete tests are better than no tests
bool testSegment() {
    // General arrangement
    const Ray<Space::Object> ray{ {1.f, 1.f}, { 2.f, 1.f } };
    const Segment seg{ {2.f, 1.f}, {2.f, 3.f} };
    const SegmentHit res = seg.intersect(ray);

    // Collinear
    const Ray<Space::Object> ray1{ {1.f,0.f}, {2.f,1.f} };
    const Segment seg1{ {5.f,2.f}, {9.f,4.f} };
    const SegmentHit res1 = seg1.intersect(ray1);

    // Parallel, but not collinear
    const Ray<Space::Object> ray2{ {1.f,10.f}, {2.f,1.f} };
    const SegmentHit res2 = seg1.intersect(ray2);
    return (res && approx(res.u, 0.25f) && !res1 && !res2);
}

bool testBoxIntersect() {
    Box<Space::Object> box;
    box.expand({ 1.f, 1.f });
    box.expand({ 4.f, 5.f });

    // Horizontal, hit
    const float t = box.intersect({ {5.f, 3.f}, {-1.f,0.f} });
    // Horizontal, inside
    const float t1 = box.intersect({ {3.f, 2.f}, {-1.f, 0.f} });
    // Vertical, miss
    const float t2 = box.intersect({ {0.f, 4.f}, {-1.f, 0.f} });
    // Random, hit
    const float t3 = box.intersect({ {0.f, 0.f}, {1.f, 2.f} });

    return (approx(t, 1.f) && approx(t1, 2.f) && (t2 == inf) && approx(t3, std::sqrt(5.f)));
}

bool testTransform() {
    const float phi = 0.25f*3.1415926f;
    const float sinPhi = std::sin(phi);
    const float cosPhi = std::cos(phi);
    // Some rotation, followed by translation
    const Transform<Space::Object, Space::World> tm1{ {{cosPhi, -sinPhi}, {sinPhi, cosPhi}}, {-3.f, 5.f} };
    const Transform<Space::World, Space::Object> itm1 = inverse(tm1);
    for (int i = 0; i < 4; ++i) {
        const Point<Space::Object> p0{ ((i == 0 || i == 3) ? 1.f : -1.f), (i < 2 ? 1.f : -1.f) };
        const Point<Space::Object> p1 = itm1 * (tm1 * p0);
        if (!approx(p0.x, p1.x) || !approx(p0.y, p1.y)) {
            return false;
        }
    }
    return true;
}

bool testTriangle() {
    // Construct a "random" triangle, point a ray at a specific point & check whether the results match
    const Ray3D::Point v0{ 1.f, -1.f, 0.f };
    const Ray3D::Point v1{ 2.f, 2.f, 1.f };
    const Ray3D::Point v2{ -3.f, 1.f, -2.f };
    const Ray3D::Vector e1 = v1 - v0;
    const Ray3D::Vector e2 = v2 - v0;
    const float origU = 0.25f, origV = 0.5f;
    const Ray3D::Point hit = v0 + origU*e1 + origV*e2;
    const Ray3D::Vector dir = { 1.f, 2.f, -1.f };
    const float dirLen = std::sqrt(dot(dir, dir));
    const Ray3D::Ray ray{ hit - dir, dir };
    const Ray3D::TriangleHit res = Ray3D::Triangle{ v0, e1, e2 }.intersect(ray);
    return (res && approx(res.t, dirLen) && approx(res.u, origU) && approx(res.v, origV));
}

int main() {
    vassert(testSegment());
    vassert(testBoxIntersect());
    vassert(testTransform());
    vassert(testTriangle());
    std::cout << "Hello, world!\n";
}