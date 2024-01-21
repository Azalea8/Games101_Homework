//
// Created by goksu on 2/25/20.
//

#include <fstream>
#include "Scene.hpp"
#include "Renderer.hpp"
#include<thread>
#include<mutex>

inline float deg2rad(const float& deg) { return deg * M_PI / 180.0; }

const float EPSILON = 0.00001;

// The main render function. This where we iterate over all pixels in the image,
// generate primary rays and cast these rays into the scene. The content of the
// framebuffer is saved to a file.
void Renderer::Render(const Scene& scene)
{
    std::vector<Vector3f> framebuffer(scene.width * scene.height);

    float scale = tan(deg2rad(scene.fov * 0.5));
    float imageAspectRatio = scene.width / (float)scene.height;

    Vector3f eye_pos(278, 273, -800);
    int m = 0;

    // 之前的实现，单线程
/*    // change the spp value to change sample ammount
    int spp = 16;
    int l = 5;
    std::cout << "SPP: " << spp << "\n";
    for (uint32_t j = 0; j < scene.height; ++j) {
        for (uint32_t i = 0; i < scene.width; ++i) {

            for(int i2 = 1; i2 < l; i2++){
                for(int j2 = 1; j2 < l; j2++){
                    // generate primary ray direction
                    float x = (2 * (i + (1.0f * i2 / l)) / (float)scene.width - 1) *
                            imageAspectRatio * scale;
                    float y = (1 - 2 * (j + (1.0f * j2 / l)) / (float)scene.height) * scale;

                    Vector3f dir = normalize(Vector3f(-x, y, 1));
                    framebuffer[m] += scene.castRay(Ray(eye_pos, dir), 0) / spp;
                }
            }

            m++;
        }
        UpdateProgress(j / (float)scene.height);
    }
    UpdateProgress(1.f);*/

    // change the spp value to change sample ammount
    int spp = 16;
    int l = 5;

    int thread_num = 16; //我的电脑有 16个逻辑处理器，所以开 16个线程。注：屏幕的高度一定要是线程数的倍数

    int thread_height = scene.height / thread_num;
    std::vector<std::thread> threads(thread_num);
    std::cout << "SPP: " << spp << "\n";

    //多线程实现
    std::mutex mtx;
    float process=0;
    float Reciprocal_Scene_height=1.f / (float)scene.height;
    auto castRay = [&](int thread_index)
    {
        int height = thread_height * (thread_index + 1);
        for (uint32_t j = height - thread_height; j < height; j++)
        {
            for (uint32_t i = 0; i < scene.width; ++i) {
                
                int m = j*scene.width+i;
                // 更改了框架
                // 对像素多次采样平均的时候，采样的都是像素中心，降噪完全依靠后续的随机数。
                // 并且超采样可以抗锯齿，但是会导致图像有点模糊，只能提高真实采样，增大图片的像素
                // 闫老师的 pdf那一页就可以明显看到，光线是从同一个像素的不同点出去的，ssp=81的时候就可以明显接近真实照片
                for(int i2 = 1; i2 < l; i2++){
                    for(int j2 = 1; j2 < l; j2++){
                        // generate primary ray direction
                        float x = (2 * (i + (1.0f * i2 / l)) / (float)scene.width - 1) *
                                imageAspectRatio * scale;
                        float y = (1 - 2 * (j + (1.0f * j2 / l)) / (float)scene.height) * scale;

                        Vector3f dir = normalize(Vector3f(-x, y, 1));

                        framebuffer[m] += scene.castRay(Ray(eye_pos, dir), 0) / (float)spp;

                        // for (int k = 0; k < 16;k++) {
                        //     // 同一个像素超采样可以不采用算数平均，像素中心的贡献值应该比靠近边缘的大，可以采用离散的正态分布作为权重
                        //     framebuffer[m] += scene.castRay(Ray(eye_pos, dir), 0) / (float)spp / 16.f;
                        // }
                    }
                }
            }
            mtx.lock();
            process = process + Reciprocal_Scene_height;
            UpdateProgress(process);
            mtx.unlock();
        }
    };

    for (int k = 0; k < thread_num; k++)
    {
        threads[k] = std::thread(castRay,k);
    }
    for (int k = 0; k < thread_num; k++)
    {
        threads[k].join();
    }
    UpdateProgress(1.f);

    // save framebuffer to file
    FILE* fp = fopen("binary.ppm", "wb");
    (void)fprintf(fp, "P6\n%d %d\n255\n", scene.width, scene.height);
    for (auto i = 0; i < scene.height * scene.width; ++i) {
        static unsigned char color[3];
        color[0] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].x), 0.6f));
        color[1] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].y), 0.6f));
        color[2] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].z), 0.6f));
        fwrite(color, 1, 3, fp);
    }
    fclose(fp);    
}
