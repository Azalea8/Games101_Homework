//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this -> bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this -> bvh -> Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    // 所有发光物体的面积
    float emit_area_sum = 0;
    // 遍历所有物体
    for (uint32_t k = 0; k < objects.size(); ++k) {
        // 如果这个物体是发光体
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }

    // 中值定理
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            // 总有一个发光体累加后能超过设定的 p，近似于随机一个发光体，但其实这种效果其实并不好，不知道为啥要这样写
            // 感觉其实可以对发光的面元编个号，对编号随机
            if (p <= emit_area_sum){
                // 对该光源采样
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// 这里已经没有深度了，光线会尽可能一直弹射下去
Vector3f Scene::shade(Intersection& hit_obj, Vector3f wo) const
{
    // 该束光线正中光源，无需着色，直接返回
    if (hit_obj.m->hasEmission())
    {
        return hit_obj.m->getEmission();
    }

    // 由于计算机对于浮点数的相等判断过于困难，只能设定一个偏差值
    const float epsilon = 0.0005f;

    switch (hit_obj.m -> getType()) {
        case DIFFUSE:{
            // 直接光照贡献
            // 直接光照只能来自场景中的光源，物体是一定能跟光源连上线的，但是其中可能有遮挡。
            // 全球面去赌的话，开销太大了
            // 所以才让直接直接光照的光路从光源的投影那块面积出发，是否遮挡后续再处理
            Vector3f Lo_dir;
            {
                // 这里是对光源采样
                float light_pdf;
                Intersection hit_light;
                sampleLight(hit_light, light_pdf);

                Vector3f obj2Light = hit_light.coords - hit_obj.coords;     // 着色点到随机点光源的向量
                Vector3f obj2LightDir = obj2Light.normalized();     // 从着色点指向点光源

                // 检查光线是否被遮挡
                auto t = intersect(Ray(hit_obj.coords, obj2LightDir));
                if (t.distance - obj2Light.norm() > -epsilon)
                {
                    Vector3f f_r = hit_obj.m -> eval(obj2LightDir, wo, hit_obj.normal);
                    float r2 = dotProduct(obj2Light, obj2Light);
                    float cosA = std::max(.0f, dotProduct(hit_obj.normal,obj2LightDir));
                    float cosB = std::max(.0f, dotProduct(hit_light.normal,-obj2LightDir));
                    Lo_dir = hit_light.emit * f_r * cosA * cosB / r2 / light_pdf;
                }
            }

            // 间接光照贡献
            // 关于为何不可以把其他物体也当做光源处理，主要原因是这样做会使问题变得更加复杂。(把所有物体投影到着色点那个球面上，间接光照也来自于这些区域)
            // 当我们考虑其他物体作为光源时，我们需要计算这些物体对整个场景的贡献，这通常涉及到更复杂的数学计算和采样方法。
            // 而单独对光源进行面积积分可以使问题得到简化，因为我们可以更容易地控制和优化采样过程。
            Vector3f Lo_indir;
            {
                // 用随机数来确保递归停止，递归层数越深，越接近数学期望
                if (get_random_float() < RussianRoulette)
                {
                    Vector3f dir2NextObj = hit_obj.m -> sample(wo, hit_obj.normal).normalized();

                    float pdf = hit_obj.m -> pdf(wo, dir2NextObj, hit_obj.normal);
                    if (pdf > epsilon)
                    {
                        Intersection nextObj = intersect(Ray(hit_obj.coords, dir2NextObj));
                        // 保证碰撞发生，并且物体不是发光体，间接光照递归
                        if (nextObj.happened && !nextObj.m->hasEmission())
                        {
                            Vector3f f_r = hit_obj.m -> eval(dir2NextObj, wo, hit_obj.normal); //BRDF
                            float cos = std::max(.0f, dotProduct(dir2NextObj, hit_obj.normal));
                            // 这里要记得除以 RussianRoulette，保证数学期望正确
                            Lo_indir = shade(nextObj, -dir2NextObj) * f_r * cos / pdf / RussianRoulette;
                        }
                    }
                }
            }

            return Lo_dir + Lo_indir;
        }
        case MIRROR:{
            Vector3f Lo_indir;
            {
                // 用随机数来确保递归停止，递归层数越深，越接近数学期望
                if (get_random_float() < RussianRoulette)
                {
                    Vector3f dir2NextObj = hit_obj.m -> sample(wo, hit_obj.normal).normalized();

                    float pdf = hit_obj.m -> pdf(wo, dir2NextObj, hit_obj.normal);
                    if (pdf > epsilon)
                    {
                        Intersection nextObj = intersect(Ray(hit_obj.coords, dir2NextObj));
                        // 保证碰撞发生
                        if(nextObj.happened) {
                            Vector3f f_r = hit_obj.m -> eval(dir2NextObj, wo, hit_obj.normal); //BRDF
                            float cos = std::max(.0f, dotProduct(dir2NextObj, hit_obj.normal));

                            // 如果不是发光体
                            if(!nextObj.m->hasEmission()) {
                                // 这里要记得除以 RussianRoulette，保证数学期望正确
                                Lo_indir = shade(nextObj, -dir2NextObj) * f_r * cos / pdf / RussianRoulette;
                            }else {
                                // TODO: 镜面反射后进入光源，按理说会是高光，不知道怎么做高光
                                // 如果不分支的话 Lo_indir = nextObj.m->getEmission() * f_r * cos / pdf / RussianRoulette;
                                // 但这里我想的是属于直接光照，没必要去赌了
                                Lo_indir = nextObj.m -> getEmission() * f_r * cos;
                            }
                        }
                    }
                }
            }

            return  Lo_indir;
        }
    }
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here
    auto hitObj = intersect(ray);
    if (!hitObj.happened) return {};
    return shade(hitObj,-ray.direction);
}