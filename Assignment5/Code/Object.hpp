#pragma once

#include "Vector.hpp"
#include "global.hpp"

class Object
{
public:
    // Sphere 和 MeshTriangle的父类
    Object()
        : materialType(DIFFUSE_AND_GLOSSY)
        , ior(1.6)
        , Kd(0.8)
        , Ks(0.2)
        , diffuseColor(0.2)
        , specularExponent(25)
    {}

    virtual ~Object() = default;

    virtual bool intersect(const Vector3f&, const Vector3f&, float&, uint32_t&, Vector2f&) const = 0;

    virtual void getSurfaceProperties(const Vector3f&, const Vector3f&, const uint32_t&, const Vector2f&, Vector3f&,
                                      Vector2f&) const = 0;

    virtual Vector3f evalDiffuseColor(const Vector2f&) const
    {
        return diffuseColor;
    }

    // material properties
    MaterialType materialType;
    float ior; // 折射率
    float Kd, Ks;
    Vector3f diffuseColor;
    float specularExponent;
};
