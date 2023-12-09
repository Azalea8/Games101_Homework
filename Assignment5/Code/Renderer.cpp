#include <fstream>
#include "Vector.hpp"
#include "Renderer.hpp"
#include "Scene.hpp"
#include <optional>

// 角度转弧度
inline float deg2rad(const float &deg){
    return deg * M_PI/180.0;
}

// 计算折射光
Vector3f reflect(const Vector3f &I, const Vector3f &N)
{
    return I - 2 * dotProduct(I, N) * N;
}


// 计算折射光
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

// 菲涅尔反射系数
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


// 如果射线与物体相交，则返回 true，否则返回 false。
//
// orig是射线的起点
// dir是射线的方向
// objects是场景中包含的物体列表
std::optional<hit_payload> trace(
        const Vector3f &orig, const Vector3f &dir,
        const std::vector<std::unique_ptr<Object>> &objects)
{
    // tNear初始为最大
    float tNear = kInfinity;

    // 碰撞信息
    std::optional<hit_payload> payload;

    // 遍历所有在场景中的物体
    for (const auto & object : objects)
    {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        // 条件 tNearK < tNear 找到离观测点最近的物体
        // intersect函数是判断光线是否与该物体相交，不同的物体函数实现不一样
        if (object -> intersect(orig, dir, tNearK, indexK, uvK) && tNearK < tNear)
        {
            // 需要清除之前写入过的信息
            payload.emplace();
            // 写入碰撞信息
            payload->hit_obj = object.get();
            payload->tNear = tNearK;
            payload->index = indexK;
            payload->uv = uvK;
            // 更新新的 tNear，确保找到离观测点最近的物体
            tNear = tNearK;
        }
    }

    // 此时返回的必定是距离观测点最近的物体的碰撞信息
    // 如果压根没有物体与光线碰撞，则返回初始的值
    return payload;
}


// 实现 Whitted风格的光传输算法（E [S*] (D|G) L）
// 这个函数计算一个由位置和方向定义的射线在交点处的颜色。注意，这个函数是递归的（它调用自身）。
//
// 如果交点处物体的材质是反射性的或者既有反射性又有折射性，那么我们计算 反射/折射 方向，并通过递归调用 castRay()函数在场景中投射两条新的射线。
// 当表面是透明的时候，我们使用菲涅尔方程的结果（它根据表面法线、入射视线方向和表面折射率计算反射和折射的比例）来混合反射和折射的颜色。
//
// 如果表面是 漫反射/高光 的，我们使用 Phong光照模型来计算交点处的颜色。
Vector3f castRay(const Vector3f &orig, const Vector3f &dir, const Scene& scene, int depth)
{
    // 超过递归深度
    if (depth > scene.maxDepth) {
        return Vector3f(0.0,0.0,0.0);
    }

    // 默认为背景颜色
    Vector3f hitColor = scene.backgroundColor;

    // 如果光线撞击到了之前我们加入到场景中的物体，则进入分支
    // 否则该像素块为背景
    if (auto payload = trace(orig, dir, scene.get_objects()); payload)
    {
        // 有了之前通过 trace函数得到的碰撞信息，算出碰撞点
        Vector3f hitPoint = orig + dir * payload->tNear;

        Vector3f N; // 法线
        Vector2f st; // 纹理坐标

        // payload中保存了碰撞的是那个物体
        // getSurfaceProperties用来求得碰撞点在物体表面的其他信息，不同物体实现不同
        payload -> hit_obj -> getSurfaceProperties(hitPoint, dir, payload->index, payload->uv, N, st);

        // 碰撞物体材质不同
        switch (payload -> hit_obj -> materialType) {
            // 既能折射又能反射
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

                // 递归调用
                Vector3f reflectionColor = castRay(reflectionRay_Orig, reflection_Direction, scene, depth + 1);
                Vector3f refractionColor = castRay(refractionRay_Orig, refraction_Direction, scene, depth + 1);

                // 使用了菲涅尔反射系数 kr来加权反射颜色和折射颜色
                float kr = fresnel(dir, N, payload -> hit_obj -> ior);
                hitColor = reflectionColor * kr + refractionColor * (1 - kr);
                break;
            }
            // 反射
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
            // 在默认情况下，我们使用 Phong光照模型。Phong模型由漫反射和镜面反射两个组成部分组成。
            default:
            {
                Vector3f lightAmt = 0, specularColor = 0;
                // 判断碰撞点是否处于阴影下，还需要从碰撞点射向光源
                Vector3f shadowPointOrig = (dotProduct(dir, N) < 0) ?
                                           hitPoint + N * scene.epsilon :
                                           hitPoint - N * scene.epsilon;

                // 遍历场景中的所有光源，并将它们的贡献相加
                // 我们还应用 Lambert余弦定律
                for (auto& light : scene.get_lights()) {
                    // 目前 lightDir为光源到碰撞点的距离
                    Vector3f lightDir = light -> position - hitPoint;
                    float light_Distance_2 = dotProduct(lightDir, lightDir);

                    // 现在 lightDir为到光源的单位向量
                    lightDir = normalize(lightDir);

                    // Phong光照模型的漫反射相关系数
                    float LdotN = std::max(0.f, dotProduct(lightDir, N));

                    // 从该点射向光源
                    // (如果碰撞到了场景中的其他物体) && (该点与碰撞物体距离 比 该点到光源的距离更近)，则说明该点处于阴影处
                    auto shadow_res = trace(shadowPointOrig, lightDir, scene.get_objects());
                    bool inShadow = shadow_res && (shadow_res -> tNear * shadow_res -> tNear < light_Distance_2);

                    // 这个地方阴影有些硬了，没有考虑全局光照
                    lightAmt += inShadow ? 0 : light -> intensity * LdotN;

                    // 反射，注意 lightDir调转方向 (原来是从碰撞点指向光源 )
                    Vector3f reflection_Direction = reflect(-lightDir, N);

                    // Phong光照模型，镜面反射的颜色贡献
                    specularColor += powf(std::max(0.f, -dotProduct(reflection_Direction, dir)),
                        payload->hit_obj->specularExponent) * light->intensity;
                }

                // payload->hit_obj->Kd, payload->hit_obj->Ks 用来线性插值
                // 前者是漫反射，后者是镜面反射
                hitColor = (lightAmt * payload->hit_obj->evalDiffuseColor(st) * payload->hit_obj->Kd)
                        + (specularColor * payload->hit_obj->Ks);
                break;
            }
        }
    }

    return hitColor;
}


// 主要的渲染函数。
// 在这里，我们遍历图像中的所有像素，生成主要光线并将这些光线投射到场景中。
// 帧缓冲区的内容被保存到文件中。
void Renderer::Render(const Scene& scene)
{
    // 创建帧缓存
    std::vector<Vector3f> framebuffer(scene.width * scene.height);


    float scale = std::tan(deg2rad(scene.fov * 0.5f));
    float imageAspectRatio = scene.width / (float)scene.height;

    // Use this variable as the eye position to start your rays.
    // eye_pos = [0, 0, 0]
    Vector3f eye_pos(0);
    int m = 0;
    // 这里更改了框架，个人习惯 i遍历行，j遍历列
    for (int i = 0; i < scene.height; ++i)
    {
        for (int j = 0; j < scene.width; ++j)
        {
            // 坐标系转换
            float x = (2 * ((j + 0.5f) / scene.width) - 1.0f) * imageAspectRatio * 1 * scale;
            float y = (1.0f - 2 * ((i + 0.5f) / scene.height)) * 1 * scale;

            // Vector3f(x, y, -1) 说明 scene在 z = -1，故 z-near距离人眼距离为 1
            Vector3f dir = normalize(Vector3f(x, y, -1));   // Don't forget to normalize this direction!

            // 这里光线追踪返回该像素点对应的颜色，记录到帧缓存中
            // 该作业的光线追踪并不是全局光照，光线如果撞击到物体就会立刻返回颜色
            // 如果是反射和折射才会递归下去
            framebuffer[m++] = castRay(eye_pos, dir, scene, 0);
        }
        // 这里是用来显示进度条，与光线追踪无关无关
        UpdateProgress(i / (float)scene.height);
    }

    // save framebuffer to file
    // 放弃使用 opencv库，反正我们只需维护 scene.height * scene.width大小的 Vector3f数组，
    // 不考虑实时显示的话，opencv对于我们来说太鸡肋
    // .ppm后缀的图像，windows下无法直接打开，vscode有个插件可以查看,不过如果图片过大，无法加载。
    // Linux下可以直接查看
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

    std::cout << std::endl << "over" << std::endl;
}
