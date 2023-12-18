#include "Renderer.hpp"
#include "Scene.hpp"
#include "Triangle.hpp"
#include "Sphere.hpp"
#include "Vector.hpp"
#include "global.hpp"
#include <chrono>

// In the main function of the program, we create the scene (create objects and
// lights) as well as set the options for the render (image width and height,
// maximum recursion depth, field-of-view, etc.). We then call the render
// function().
int main(int argc, char** argv)
{

    // Change the definition here to change resolution
    Scene scene(784, 784);

    Material* red = new Material(DIFFUSE, Vector3f(0.0f));
    red->Kd = Vector3f(0.63f, 0.065f, 0.05f);

    Material* green = new Material(DIFFUSE, Vector3f(0.0f));
    green->Kd = Vector3f(0.14f, 0.45f, 0.091f);

    Material* white = new Material(DIFFUSE, Vector3f(0.0f));
    white->Kd = Vector3f(0.725f, 0.71f, 0.68f);

    Material* white_mirror = new Material(MIRROR, Vector3f(0.0f));
    white_mirror->Kd = Vector3f(0.0f, 0.0f, 0.0f);
    white_mirror->ior = 40.0f;

    Material* microfacet = new Material(MICROFACET, Vector3f(0.0f));
    microfacet->Kd = Vector3f(0.3, 0.3, 0.25);
    microfacet->roughness = 0.5f;
    // microfacet->ior = 12.0f;

    Material* micro_Al = new Material(MICROFACET, Vector3f(0.0f));
    micro_Al->Kd = Vector3f(0.913f, 0.921f, 0.925f);
    micro_Al->roughness = 0.001f;
    micro_Al->metalness = 1.0f;

    Material* Au = new Material(DIFFUSE, Vector3f(0.0f));
    Au->Kd = Vector3f(1.000f, 0.766f, 0.336f);

    Material* micro_Au = new Material(MICROFACET, Vector3f(0.0f));
    micro_Au->Kd = Vector3f(1.000f, 0.766f, 0.336f);
    micro_Au->roughness = 0.5f;
    micro_Au->metalness = 0.8f;

    Material* light = new Material(DIFFUSE, (8.0f * Vector3f(0.747f+0.058f, 0.747f+0.258f, 0.747f) + 15.6f * Vector3f(0.740f+0.287f,0.740f+0.160f,0.740f) + 18.4f *Vector3f(0.737f+0.642f,0.737f+0.159f,0.737f)));
    light->Kd = Vector3f(0.65f);

    MeshTriangle floor("../models/cornellbox/floor.obj", white);
    MeshTriangle shortbox("../models/cornellbox/shortbox.obj", micro_Au);
    MeshTriangle tallbox("../models/cornellbox/tallbox.obj", micro_Al);
    MeshTriangle left("../models/cornellbox/left.obj", red);
    MeshTriangle right("../models/cornellbox/right.obj", green);
    MeshTriangle light_("../models/cornellbox/light.obj", light);
    // MeshTriangle bunny("../models/bunny/bunny.obj", white_mirror, Vector3f(200,-60,150), Vector3f(1500,1500,1500), Vector3f(-1,0,0), Vector3f(0,1,0), Vector3f(0,0,-1));
    // Sphere sphere(Vector3f(150, 100, 200), 100, micro_Au);

    scene.Add(&floor);
    scene.Add(&shortbox);
    scene.Add(&tallbox);
    scene.Add(&left);
    scene.Add(&right);
    scene.Add(&light_);
    // scene.Add(&bunny);
    // scene.Add(&sphere);
    
    scene.buildBVH();

    Renderer r;

    auto start = std::chrono::system_clock::now();
    r.Render(scene);
    auto stop = std::chrono::system_clock::now();

    std::cout << "Render complete: \n";
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::hours>(stop - start).count() << " hours\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::minutes>(stop - start).count() << " minutes\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds\n";

    return 0;
}