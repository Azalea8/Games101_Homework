#pragma once
#include "Scene.hpp"

// 碰撞信息
// tNear存储到最近相交物体的距离。
// index存储相交物体为网格时的相交三角形的索引。
// uv存储相交点的 u和 v重心坐标。
// *hit_Obj存储指向相交物体的指针（用于获取材质信息等）。
struct hit_payload
{
    float tNear;
    uint32_t index;
    Vector2f uv;
    Object* hit_obj;
};

class Renderer
{
public:
    void Render(const Scene& scene);

private:
};