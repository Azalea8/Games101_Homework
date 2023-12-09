#include "Scene.hpp"
#include "Sphere.hpp"
#include "Triangle.hpp"
#include "Light.hpp"
#include "Renderer.hpp"

// In the main function of the program, we create the scene (create objects and lights)
// as well as set the options for the render (image width and height, maximum recursion
// depth, field-of-view, etc.). We then call the render function().
int main()
{
    Scene scene(2000, 2000);

    // 创建一个球，该球继承 Object类
    auto sph1 = std::make_unique<Sphere>(Vector3f(-1, 0, -12), 2);
    // 球的材质
    sph1 -> materialType = DIFFUSE_AND_GLOSSY;
    // 真实的颜色
    sph1 -> diffuseColor = Vector3f(1.0, 1.0, 1.0);

    auto sph2 = std::make_unique<Sphere>(Vector3f(0.5, -0.5, -8), 1.5);
    // 球的折射率
    sph2 -> ior = 1.5;
    sph2 -> materialType = REFLECTION_AND_REFRACTION;

    // 将物体加入到场景中
    scene.Add(std::move(sph1));
    scene.Add(std::move(sph2));

    // 定义了两个三角形，由于有公用顶点的存在
    // vertIndex来解释哪些顶点组成一个三角形，此处为 [0, 1, 3] 和 [1, 2, 3]
    Vector3f verts[4] = {{-5,-3,-6}, {5,-3,-6}, {5,-3,-16}, {-5,-3,-16}};
    uint32_t vertIndex[6] = {0, 1, 3, 1, 2, 3};

    // 纹理坐标
    Vector2f st[4] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};

    // 创建一堆三角形
    auto mesh = std::make_unique<MeshTriangle>(verts, vertIndex, 2, st);
    mesh -> materialType = DIFFUSE_AND_GLOSSY;

    scene.Add(std::move(mesh));

    scene.Add(std::make_unique<Light>(Vector3f(-20, 70, 20), 0.5));
    scene.Add(std::make_unique<Light>(Vector3f(30, 50, -12), 0.5));    

    // 准备工作就绪，渲染开始
    Renderer r;
    r.Render(scene);

    return 0;
}