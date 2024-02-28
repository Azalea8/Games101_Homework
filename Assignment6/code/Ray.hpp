//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_RAY_H
#define RAYTRACING_RAY_H
#include "Vector.hpp"
struct Ray{
    //Destination = origin + t * direction
    Vector3f origin; // 起点
    Vector3f direction; // 方向向量
    Vector3f direction_inv; // 倒数

    double t;   //transportation time,
    double t_min, t_max;

    Ray(const Vector3f& ori, const Vector3f& dir, const double _t = 0.0): origin(ori), direction(dir), t(_t) {
        direction_inv = Vector3f(1.f/direction.x, 1.f/direction.y, 1.f/direction.z);
        t_min = 0.0;
        t_max = std::numeric_limits<double>::max();

    }

    Vector3f operator()(double t) const{return origin + direction * t;}

    friend std::ostream &operator<<(std::ostream& os, const Ray& r){
        os<<"[origin:="<<r.origin<<", direction="<<r.direction<<", time="<< r.t<<"]\n";
        return os;
    }
};
#endif //RAYTRACING_RAY_H
