#include <fstream>
#include "Vector.hpp"
#include "Renderer.hpp"
#include "Scene.hpp"
#include <optional>

// 内联函数,减少开销
// https://zhuanlan.zhihu.com/p/151995167
inline float deg2rad(const float &deg){
    return deg * M_PI/180.0;
}

// Compute reflection direction
Vector3f reflect(const Vector3f &I, const Vector3f &N)
{
    return I - 2 * dotProduct(I, N) * N;
}

// [comment]
// Compute refraction direction using Snell's law
//
// We need to handle with care the two possible situations:
//
//    - When the ray is inside the object
//
//    - When the ray is outside.
//
// If the ray is outside, you need to make cosi positive cosi = -N.I
//
// If the ray is inside, you need to invert the refractive indices and negate the normal N
// [/comment]
Vector3f refract(const Vector3f &I, const Vector3f &N, const float &ior)
{
    float cosi = clamp(-1, 1, dotProduct(I, N));
    float etai = 1, etat = ior;
    Vector3f n = N;
    if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -N; }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
}

// [comment]
// Compute Fresnel equation
//
// \param I is the incident view direction
//
// \param N is the normal at the intersection point
//
// \param ior is the material refractive index
// [/comment]
float fresnel(const Vector3f &I, const Vector3f &N, const float &ior)
{
    float cosi = clamp(-1, 1, dotProduct(I, N));
    float etai = 1, etat = ior;
    if (cosi > 0) {  std::swap(etai, etat); }
    // Compute sini using Snell's law
    float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
    // Total internal reflection
    if (sint >= 1) {
        return 1;
    }
    else {
        float cost = sqrtf(std::max(0.f, 1 - sint * sint));
        cosi = fabsf(cosi);
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        return (Rs * Rs + Rp * Rp) / 2;
    }
    // As a consequence of the conservation of energy, transmittance is given by:
    // kt = 1 - kr;
}

// [comment]
// 如果射线与物体相交，则返回 true，否则返回 false。
//
// \param orig是射线的起点
// \param dir是射线的方向
// \param objects是场景中包含的物体列表
// \param[out] tNear存储到最近相交物体的距离。
// \param[out] index存储相交物体为网格时的相交三角形的索引。
// \param[out] uv存储相交点的 u和v重心坐标。
// \param[out] *hitObject存储指向相交物体的指针（用于获取材质信息等）。
// \param isShadowRay是否为阴影射线。如果我们找到了一个相交点，我们可以更早地从函数中返回。
// [/comment]
std::optional<hit_payload> trace(
        const Vector3f &orig, const Vector3f &dir,
        const std::vector<std::unique_ptr<Object>> &objects)
{
    float tNear = kInfinity;
    std::optional<hit_payload> payload;
    for (const auto & object : objects)
    {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (object -> intersect(orig, dir, tNearK, indexK, uvK) && tNearK < tNear)
        {
            payload.emplace();
            payload->hit_obj = object.get();
            payload->tNear = tNearK;
            payload->index = indexK;
            payload->uv = uvK;
            tNear = tNearK;
        }
    }

    return payload;
}

// [comment]
// 实现 Whitted风格的光传输算法（E [S*] (D|G) L）
//
// 这个函数计算一个由位置和方向定义的射线在交点处的颜色。注意，这个函数是递归的（它调用自身）。
//
// 如果交点处物体的材质是反射性的或者既有反射性又有折射性，那么我们计算反射/折射方向，并通过递归调用castRay()函数在场景中投射两条新的射线。
// 当表面是透明的时候，我们使用菲涅尔方程的结果（它根据表面法线、入射视线方向和表面折射率计算反射和折射的比例）来混合反射和折射的颜色。
//
// 如果表面是漫反射/高光的，我们使用Phong光照模型来计算交点处的颜色。
// [/comment]
Vector3f castRay(const Vector3f &orig, const Vector3f &dir, const Scene& scene, int depth)
{
    if (depth > scene.maxDepth) {
        return Vector3f(0.0,0.0,0.0);
    }

    Vector3f hitColor = scene.backgroundColor;

    // 如果光线撞击到了之前我们加入到场景中的物体，则进入分支
    if (auto payload = trace(orig, dir, scene.get_objects()); payload)
    {
        Vector3f hitPoint = orig + dir * payload->tNear;
        Vector3f N; // normal
        Vector2f st; // st coordinates
        payload -> hit_obj -> getSurfaceProperties(hitPoint, dir, payload->index, payload->uv, N, st);
        switch (payload -> hit_obj -> materialType) {
            case REFLECTION_AND_REFRACTION:
            {
                // 反射光
                Vector3f reflection_Direction = normalize(reflect(dir, N));

                // 折射光
                Vector3f refraction_Direction = normalize(refract(dir, N, payload -> hit_obj->ior));

                // 判断光线是从内部射向表面，还是外部射向表面
                // 这个偏移量的目的是避免光线与命中点所在的表面相交
                Vector3f reflectionRay_Orig = (dotProduct(reflection_Direction, N) < 0) ?
                                             hitPoint - N * scene.epsilon :
                                             hitPoint + N * scene.epsilon;
                Vector3f refractionRay_Orig = (dotProduct(refraction_Direction, N) < 0) ?
                                             hitPoint - N * scene.epsilon :
                                             hitPoint + N * scene.epsilon;
                Vector3f reflectionColor = castRay(reflectionRay_Orig, reflection_Direction, scene, depth + 1);
                Vector3f refractionColor = castRay(refractionRay_Orig, refraction_Direction, scene, depth + 1);
                float kr = fresnel(dir, N, payload -> hit_obj -> ior);
                hitColor = reflectionColor * kr + refractionColor * (1 - kr);
                break;
            }
            case REFLECTION:
            {
                float kr = fresnel(dir, N, payload->hit_obj->ior);
                Vector3f reflectionDirection = reflect(dir, N);
                Vector3f reflectionRayOrig = (dotProduct(reflectionDirection, N) < 0) ?
                                             hitPoint + N * scene.epsilon :
                                             hitPoint - N * scene.epsilon;
                hitColor = castRay(reflectionRayOrig, reflectionDirection, scene, depth + 1) * kr;
                break;
            }
            default:
            {
                // [comment]
                // We use the Phong illumation model int the default case. The phong model
                // is composed of a diffuse and a specular reflection component.
                // [/comment]
                Vector3f lightAmt = 0, specularColor = 0;
                Vector3f shadowPointOrig = (dotProduct(dir, N) < 0) ?
                                           hitPoint + N * scene.epsilon :
                                           hitPoint - N * scene.epsilon;
                // [comment]
                // Loop over all lights in the scene and sum their contribution up
                // We also apply the lambert cosine law
                // [/comment]
                for (auto& light : scene.get_lights()) {
                    Vector3f lightDir = light->position - hitPoint;
                    // square of the distance between hitPoint and the light
                    float lightDistance2 = dotProduct(lightDir, lightDir);
                    lightDir = normalize(lightDir);
                    float LdotN = std::max(0.f, dotProduct(lightDir, N));
                    // is the point in shadow, and is the nearest occluding object closer to the object than the light itself?
                    auto shadow_res = trace(shadowPointOrig, lightDir, scene.get_objects());
                    bool inShadow = shadow_res && (shadow_res->tNear * shadow_res->tNear < lightDistance2);

                    lightAmt += inShadow ? 0 : light->intensity * LdotN;
                    Vector3f reflection_Direction = reflect(-lightDir, N);

                    specularColor += powf(std::max(0.f, -dotProduct(reflection_Direction, dir)),
                        payload->hit_obj->specularExponent) * light->intensity;
                }

                hitColor = lightAmt * payload->hit_obj->evalDiffuseColor(st) * payload->hit_obj->Kd + specularColor * payload->hit_obj->Ks;
                break;
            }
        }
    }

    return hitColor;
}

// [comment]
// The main render function. This where we iterate over all pixels in the image, generate
// primary rays and cast these rays into the scene. The content of the framebuffer is
// saved to a file.
// [/comment]
void Renderer::Render(const Scene& scene)
{
    std::vector<Vector3f> framebuffer(scene.width * scene.height);

    float scale = std::tan(deg2rad(scene.fov * 0.5f));
    float imageAspectRatio = scene.width / (float)scene.height;

    // Use this variable as the eye position to start your rays.
    Vector3f eye_pos(0);
    int m = 0;
    // 这里更改了框架，个人习惯 i遍历行，j遍历列
    for (int i = 0; i < scene.height; ++i)
    {
        for (int j = 0; j < scene.width; ++j)
        {
            // generate primary ray direction
            float x = (2 * ((j + 0.5f) / scene.width) - 1.0f) * imageAspectRatio * 1 * scale;
            float y = (1.0f - 2 * ((i + 0.5f) / scene.height)) * 1 * scale;

            // Vector3f(x, y, -1) 说明 scene在 z = -1，故 z-near距离人眼距离为 1
            Vector3f dir = normalize(Vector3f(x, y, -1)); // Don't forget to normalize this direction!
            framebuffer[m++] = castRay(eye_pos, dir, scene, 0);
        }
        UpdateProgress(i / (float)scene.height);
    }

    // save framebuffer to file
    FILE* fp = fopen("binary.ppm", "wb");
    (void)fprintf(fp, "P6\n%d %d\n255\n", scene.width, scene.height);
    for (auto i = 0; i < scene.height * scene.width; ++i) {
        static unsigned char color[3];
        color[0] = (char)(255 * clamp(0, 1, framebuffer[i].x));
        color[1] = (char)(255 * clamp(0, 1, framebuffer[i].y));
        color[2] = (char)(255 * clamp(0, 1, framebuffer[i].z));
        fwrite(color, 1, 3, fp);
    }
    fclose(fp);    
}
