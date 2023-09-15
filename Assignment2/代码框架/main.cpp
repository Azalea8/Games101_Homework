// clang-format off
#include <iostream>
#include <opencv2/opencv.hpp>
#include "rasterizer.hpp"
#include "global.hpp"
#include "Triangle.hpp"

constexpr double MY_PI = 3.1415926;

// View变换
Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1,0,0,-eye_pos[0],
                 0,1,0,-eye_pos[1],
                 0,0,1,-eye_pos[2],
                 0,0,0,1;

    view = translate*view;

    return view;
}

// Model变换
Eigen::Matrix4f get_model_matrix(float angle)
{
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();
    angle = angle * MY_PI / 180.f;
    model << cos(angle), 0, sin(angle), 0,
            0, 1, 0, 0,
            -sin(angle), 0, cos(angle), 0,
            0, 0, 0, 1;
    return model;
}

// Projection变换
Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio, float zNear, float zFar)
{
    // TODO: Copy-paste your implementation from the previous assignment.
    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();
    Eigen::Matrix4f m;
    m << -zNear, 0, 0, 0,
          0, -zNear, 0, 0,
          0, 0, -zNear - zFar, -zNear * zFar,
          0, 0, 1, 0;

    float halve = (eye_fov/2) * MY_PI/180;
    float top = tan(halve) * zNear;
    float bottom = -top;
    float right = top * aspect_ratio;
    float left = -right;
    Eigen::Matrix4f n, p;
    n << 2/(right - left), 0, 0, 0,
            0, 2/(top - bottom), 0, 0,
            0, 0, 2/(zFar - zNear), 0,
            0, 0, 0, 1;

    p << 1, 0, 0, -(right + left)/2,
            0, 1, 0, -(top + bottom)/2,
            0, 0, 1, -(zFar + zNear)/2,
            0, 0, 0, 1;

    projection = n * p * m;

    return projection;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = false;
    std::string filename = "output.png";

    if (argc == 2)
    {
        command_line = true;
        filename = std::string(argv[1]);
    }

    rst::rasterizer r(700, 1400); // 光栅器初始化

    Eigen::Vector3f eye_pos = {0,0,5}; // 设置观测坐标

    // 所有顶点坐标
    std::vector<Eigen::Vector3f> pos
            {
                    {2, 0, -2},
                    {0, 2, -2},
                    {-2, 0, -2},
                    {3.5, -1, -5},
                    {2.5, 1.5, -5},
                    {-1, 0.5, -5}
            };

    // 对顶点坐标进行分组，每组为三个顶点构成一个三角形
    std::vector<Eigen::Vector3i> ind
            {
                    {0, 1, 2},
                    {3, 4, 5}
            };

    // 顶点颜色
    std::vector<Eigen::Vector3f> cols
            {
                    {217.0, 238.0, 185.0},
                    {217.0, 238.0, 185.0},
                    {217.0, 238.0, 185.0},

                    {185.0, 217.0, 238.0},
                    {185.0, 217.0, 238.0},
                    {185.0, 217.0, 238.0},
            };

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);
    auto col_id = r.load_colors(cols);

    int key = 0;
    int frame_count = 0;

    if (command_line)
    {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 2, 0.1, 50));

        r.draw(pos_id, ind_id, col_id, rst::Primitive::Triangle);

        cv::Mat image(700, 1400, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);

        cv::imwrite(filename, image);

        return 0;
    }

    while(key != 27)
    {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 2, 0.1, 50));

        r.draw(pos_id, ind_id, col_id, rst::Primitive::Triangle);

        cv::Mat image(700, 1400, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);
        cv::imshow("image", image);
        key = cv::waitKey(10);

        std::cout << "frame count: " << frame_count++ << '\n';
        // std::cout << "eye_pos_z:  " << eye_pos[2] << '\n';

        /*if (key == 'a') {
            eye_pos[0] -= 1;
        } else if (key == 'd') {
            eye_pos[0] += 1;
        } else if(key == 'w') {
            eye_pos[2] -= 1;
        } else if(key == 's') {
            eye_pos[2] += 1;
        } else if(key == 'i') {
            eye_pos[1] += 1;
        } else if(key == 'k') {
            eye_pos[1] -= 1;
        }*/

        if (key == 'a') {
            angle -= 10;
        } else if (key == 'd') {
            angle += 10;
        }
    }

    return 0;
}
// clang-format on