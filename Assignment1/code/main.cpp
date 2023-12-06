#include "Triangle.hpp"
#include "rasterizer.hpp"
#include <Eigen/Eigen>
#include <iostream>
#include <opencv2/opencv.hpp>

constexpr double MY_PI = 3.1415926;

// view变换
Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0],
                 0, 1, 0, -eye_pos[1],
                 0, 0, 1, -eye_pos[2],
                 0, 0, 0, 1;

    view = translate * view;

    return view;
}

// model变换
Eigen::Matrix4f get_model_matrix(float rotation_angle)
{
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();

    // TODO: Implement this function
    // Create the model matrix for rotating the triangle around the Z axis.
    // Then return it.

    Eigen::Matrix4f translate;
    float angle = rotation_angle * MY_PI / 180.f;
    translate << cos(angle),-sin(angle),0,0,       //绕z的旋转矩阵
                sin(angle),cos(angle),0,0,
                0,0,1,0,
                0,0,0,1;

    return translate * model;
}

// project变换
Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio,
                                      float zNear, float zFar)
{
    // Students will implement this function

    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();

    // TODO: Implement this function
    // Create the projection matrix for the given parameters.
    // Then return it.

    // 注意：传递参数时，zNear, zFar作为距离为正数，构造矩阵时，应该换算为坐标
    Eigen::Matrix4f m;
    m << -zNear, 0, 0, 0,
        0, -zNear, 0, 0,
        0, 0, -zNear - zFar, zNear * -zFar,
        0, 0, 1, 0;

    float halve = (eye_fov/2.0) * MY_PI/180;
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

int main(int argc, const char** argv) // 主函数带有参数，main函数会有分支

{
    float angle = 0;
    bool command_line = false;
    std::string filename = "output.png";

    if (argc >= 3) {
        command_line = true;
        angle = std::stof(argv[2]); // -r by default
        if (argc == 4) {
            filename = std::string(argv[3]);
        }
        else
            return 0;
    }

    // 程序入口
    rst::rasterizer r(700, 700); // 构造光栅器

    Eigen::Vector3f eye_pos = {0, 0, 5}; // 设置观测坐标，这里默认朝向-Z轴观测，不然 View变换还会涉及到旋转

    // 所有三角形的顶点坐标。这里分了组，在 opengl中无需分组，EBO一维数组一把梭，EAO负责解释哪些数据是做什么的
    std::vector<Eigen::Vector3f> pos{{2, 1, -2}, {1, 2, -2}, {-2, 1, -2}};

    std::vector<Eigen::Vector3i> ind{{0, 1, 2}, }; // 对顶点进行分组, 每组的三个值对应一个三角形

    // 将三角形信息存入光栅器，并返回唯一标识以便查询
    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);

    // std::cout << pos_id.pos_id << ' ' << ind_id.ind_id << std::endl;

    int key = 0; // 响应键盘操作
    int frame_count = 0; // 当前是第几帧

    if (command_line) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);

        cv::imwrite(filename, image);

        return 0;
    }

    while (key != 27) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth); // 清除帧缓存的图像，准备重新绘制

        r.set_model(get_model_matrix(angle)); // 向光栅器传入新的 Model变换矩阵
        r.set_view(get_view_matrix(eye_pos)); // 向光栅器传入新的View变换矩阵
        r.set_projection(get_projection_matrix(60, 1, 0.1, 50)); // 向光栅器传入新的Projection变换矩阵

        // 光栅器开始渲染。参数：所有的顶点，顶点所构成的图形顺序，图形类型。
        r.draw(pos_id, ind_id, rst::Primitive::Triangle);

        // openCV显示图像，数据存储在光栅器的帧缓存中
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::imshow("image", image);
        key = cv::waitKey(10);

        std::cout << "frame count: " << frame_count++ << '\n';
        // std::cout << "eye_pos:\n " << eye_pos << '\n';

        // 对物体在三维空间中绕 Z轴进行旋转
        if (key == 'z') {
            angle += 1;
        }
        else if (key == 'c') {
            angle -= 1;
        }

        // 观测坐标变化，不涉及观测角度变化
        if (key == 'a') {
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
        }
    }

    return 0;
}
