#pragma once

#include "Object.hpp"
#include "Vector.hpp"

class Sphere : public Object
{
public:
    // Sphere的构造函数很简陋，这里是显式曲面，判断光线是否与其相交比较简单
    Sphere(const Vector3f& c, const float& r)
        : center(c)
        , radius(r)
        , radius2(r * r)
    {}

    // 判断光线是否与球体相交
    bool intersect(const Vector3f& orig, const Vector3f& dir, float& tnear, uint32_t&, Vector2f&) const override
    {
        // 详见 Lecture-13 pdf 23页
        Vector3f L = orig - center;
        float a = dotProduct(dir, dir);
        float b = 2 * dotProduct(dir, L);
        float c = dotProduct(L, L) - radius2;
        float t0, t1;
        if (!solveQuadratic(a, b, c, t0, t1)) return false;
        if (t0 < 0)
            t0 = t1;
        if (t1 < 0)
            return false;

        // 返回最近的
        tnear = t0;

        return true;
    }

    // 碰撞点在物体表面的其他信息，如法线
    void getSurfaceProperties(const Vector3f& P, const Vector3f&, const uint32_t&, const Vector2f&,
                              Vector3f& N, Vector2f&) const override
    {
        // 显式曲面就是爽，只需要求得该点的法线就行
        // 法线从球心指向碰撞点
        N = normalize(P - center);
    }

    Vector3f center;
    float radius, radius2;
};
