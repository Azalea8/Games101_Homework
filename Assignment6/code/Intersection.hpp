//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_INTERSECTION_H
#define RAYTRACING_INTERSECTION_H
#include "Vector.hpp"
#include "Material.hpp"
class Object;
class Sphere;

struct Intersection
{
    Intersection(){
        happened=false; // 碰撞是否发生
        coords=Vector3f(); // 碰撞点在三角形内部的重心坐标
        normal=Vector3f(); // 法线
        distance= std::numeric_limits<double>::max(); // 距离，也就是 t
        obj =nullptr; // 碰撞点所属物体
        m=nullptr; // 物体的材质信息
    }
    bool happened;
    Vector3f coords;
    Vector3f normal;
    double distance;
    Object* obj;
    Material* m;
};
#endif //RAYTRACING_INTERSECTION_H
