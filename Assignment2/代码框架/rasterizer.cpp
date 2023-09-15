// clang-format off
//
// Created by goksu on 4/6/19.
//

#include <algorithm>
#include <vector>
#include "rasterizer.hpp"
#include <opencv2/opencv.hpp>
#include <cmath>

#define SSAA false

rst::pos_buf_id rst::rasterizer::load_positions(const std::vector<Eigen::Vector3f> &positions) {
    auto id = get_next_id();
    pos_buf.emplace(id, positions);

    return {id};
}

rst::ind_buf_id rst::rasterizer::load_indices(const std::vector<Eigen::Vector3i> &indices) {
    auto id = get_next_id();
    ind_buf.emplace(id, indices);

    return {id};
}

rst::col_buf_id rst::rasterizer::load_colors(const std::vector<Eigen::Vector3f> &cols) {
    auto id = get_next_id();
    col_buf.emplace(id, cols);

    return {id};
}

auto to_vec4(const Eigen::Vector3f &v3, float w = 1.0f) {
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}


static bool insideTriangle(float x, float y, const Vector3f *_v) {
    // TODO : Implement this function to check if the point (x, y) is inside the triangle represented by _v[0], _v[1], _v[2]
    Eigen::Vector3f p0p1(_v[1].x() - _v[0].x(), _v[1].y() - _v[0].y(), 1.0f);
    Eigen::Vector3f p1p2(_v[2].x() - _v[1].x(), _v[2].y() - _v[1].y(), 1.0f);
    Eigen::Vector3f p2p0(_v[0].x() - _v[2].x(), _v[0].y() - _v[2].y(), 1.0f);

    Eigen::Vector3f p0p(x - _v[0].x(), y - _v[0].y(), 1.0f);
    Eigen::Vector3f p1p(x - _v[1].x(), y - _v[1].y(), 1.0f);
    Eigen::Vector3f p2p(x - _v[2].x(), y - _v[2].y(), 1.0f);

    // std::cout << p0p1.cross(p0p) << std::endl;

    if (p0p1.cross(p0p).z() > 0.f) {
        return p1p2.cross(p1p).z() > 0.f && p2p0.cross(p2p).z() > 0.f;
    } else {
        return p1p2.cross(p1p).z() < 0.f && p2p0.cross(p2p).z() < 0.f;
    }
}

// 平滑插值, 计算权重
static std::tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector3f *v) {
    float alpha = (x * (v[1].y() - v[2].y()) + (v[2].x() - v[1].x()) * y + v[1].x() * v[2].y() - v[2].x() * v[1].y()) /
                  (v[0].x() * (v[1].y() - v[2].y()) + (v[2].x() - v[1].x()) * v[0].y() + v[1].x() * v[2].y() -
                   v[2].x() * v[1].y());
    float beta = (x * (v[2].y() - v[0].y()) + (v[0].x() - v[2].x()) * y + v[2].x() * v[0].y() - v[0].x() * v[2].y()) /
                 (v[1].x() * (v[2].y() - v[0].y()) + (v[0].x() - v[2].x()) * v[1].y() + v[2].x() * v[0].y() -
                  v[0].x() * v[2].y());
    float gamma = (x * (v[0].y() - v[1].y()) + (v[1].x() - v[0].x()) * y + v[0].x() * v[1].y() - v[1].x() * v[0].y()) /
                  (v[2].x() * (v[0].y() - v[1].y()) + (v[1].x() - v[0].x()) * v[2].y() + v[0].x() * v[1].y() -
                   v[1].x() * v[0].y());
    return {alpha, beta, gamma};
}

void rst::rasterizer::draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, col_buf_id col_buffer, Primitive type) {
    // 取出缓存在光栅器中的顶点数据
    auto &buf = pos_buf[pos_buffer.pos_id];
    auto &ind = ind_buf[ind_buffer.ind_id];
    auto &col = col_buf[col_buffer.col_id];

    float f1 = (50 - 0.1) / 2.0;
    float f2 = (50 + 0.1) / 2.0;

    Eigen::Matrix4f mvp = projection * view * model; // MVP变换

    for (auto &i: ind) {
        Triangle t;

        // v数组代表变换到屏幕的三角形，共有三个元素代表三个顶点(用齐次坐标表示)
        Eigen::Vector4f v[] = {
                mvp * to_vec4(buf[i[0]], 1.0f),
                mvp * to_vec4(buf[i[1]], 1.0f),
                mvp * to_vec4(buf[i[2]], 1.0f)
        };

        // 齐次坐标归一化使其表示一个点
        // 但是会丢失 W的值，（三维空间中的深度信息丢失）
        // 作业3 中得到修正
        for (auto &vec: v) {
            vec /= vec.w();
        }

        //视口变换
        for (auto &vert: v) {
            vert.x() = 0.5 * width * (vert.x() + 1.0);
            vert.y() = 0.5 * height * (vert.y() + 1.0);
            vert.z() = vert.z() * f1 + f2;
        }

        for (int j = 0; j < 3; ++j) {
            t.setVertex(j, v[j].head<3>());
            // t.setVertex(i, v[i].head<3>());
            // t.setVertex(i, v[i].head<3>());
        }

        auto col_x = col[i[0]];
        auto col_y = col[i[1]];
        auto col_z = col[i[2]];

        t.setColor(0, col_x[0], col_x[1], col_x[2]);
        t.setColor(1, col_y[0], col_y[1], col_y[2]);
        t.setColor(2, col_z[0], col_z[1], col_z[2]);

        rasterize_triangle(t);

        if (SSAA) {
            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    Eigen::Vector3f color(0, 0, 0);
                    for (int i = 0; i < 4; i++)
                        color += frame_buf_2xSSAA[get_index(x, y)][i];
                    color /= 4;
                    set_pixel(Eigen::Vector3f(x, y, 1.0f), color);
                }
            }
        }
    }
}

//Screen space rasterization
void rst::rasterizer::rasterize_triangle(const Triangle &t) {
    auto v = t.toVector4();

    // TODO : Find out the bounding box of current triangle.
    // iterate through the pixel and find if the current pixel is inside the triangle

    // If so, use the following code to get the interpolated z value.
    //auto[alpha, beta, gamma] = computeBarycentric2D(x, y, t.v);
    //float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
    //float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
    //z_interpolated *= w_reciprocal;

    // TODO : set the current pixel (use the set_pixel function) to the color of the triangle (use getColor function) if it should be painted.

    std::vector<float> x_arry{v[0].x(), v[1].x(), v[2].x()};
    std::vector<float> y_arry{v[0].y(), v[1].y(), v[2].y()};
    std::sort(x_arry.begin(), x_arry.end());
    std::sort(y_arry.begin(), y_arry.end());
    int x_min = floor(x_arry[0]), x_max = ceil(x_arry[2]),
            y_min = floor(y_arry[0]), y_max = ceil(y_arry[2]);

    // 光栅化
    for (int x = x_min; x < x_max; x++) {
        for (int y = y_min; y < y_max; y++) {
            // 不在视锥中的物体进行裁切
            if (x < 0 || x >= width || y < 0 || y >= height) {
                continue;
            }

            if (SSAA) {
                int index = 0;
                // 划分四个小的像素
                for (float i = 0.25; i < 1.0; i += 0.5) {
                    for (float j = 0.25; j < 1.0; j += 0.5) {
                        if (insideTriangle(x + i, y + j, t.v)) {
                            auto [alpha, beta, gamma] = computeBarycentric2D(x + i, y + j, t.v);
                            float w_reciprocal = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                            float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() +
                                                   gamma * v[2].z() / v[2].w();
                            z_interpolated *= w_reciprocal;
                            if (abs(z_interpolated) < depth_buf_2xSSAA[get_index(x, y)][index]) {
                                Eigen::Vector3f point(x + i, y + j, 1.0f);
                                frame_buf_2xSSAA[get_index(x, y)][index] = t.getColor();
                                depth_buf_2xSSAA[get_index(x, y)][index] = abs(z_interpolated);
                            }
                        }
                        index++;
                    }
                }
            } else {
                if (insideTriangle(x + 0.5f, y + 0.5f, t.v)) {
                    auto [alpha, beta, gamma] = computeBarycentric2D(x, y, t.v);

                    // 这三步为了透视矫正。 Z的插值应该是在三维空间下才准确，在二维平面下直接插值计算权重的结果不正确
                    // 透视矫正插值：https://zhuanlan.zhihu.com/p/144331875
                    // 作业1,2 框架相关问题： https://zhuanlan.zhihu.com/p/509902950
                    // w_reciprocal已经是三维空间下的正确深度信息，但对于遮挡关系还是选择二维平面下的深度信息，z_interpolated最后值为透视校正插值后的二维平面下正确坐标深度信息
                    // 当然也可以用 w-buffer，但普遍选择 z-buffer
                    float w_reciprocal = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                    float z_interpolated =
                            alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                    z_interpolated *= w_reciprocal;

                    // 注意符号！！！ 深度缓存初始化是全为正无穷，个人习惯坐标的 Z值全取负值
                    // 插值出来的 Z值与深度缓存中记录的 Z值作比较，set_pixel将帧缓存中该像素的颜色更新，最后深度缓存更新
                    if (abs(z_interpolated) < depth_buf[get_index(x, y)]) {
                        Eigen::Vector3f point(x, y, 1.0f);
                        set_pixel(point, t.getColor());
                        depth_buf[get_index(x, y)] = abs(z_interpolated);
                    }
                }
            }
        }
    }
}


void rst::rasterizer::set_model(const Eigen::Matrix4f &m) {
    model = m;
}

void rst::rasterizer::set_view(const Eigen::Matrix4f &v) {
    view = v;
}

void rst::rasterizer::set_projection(const Eigen::Matrix4f &p) {
    projection = p;
}

void rst::rasterizer::clear(rst::Buffers buff) {
    if ((buff & rst::Buffers::Color) == rst::Buffers::Color) {
        std::fill(frame_buf.begin(), frame_buf.end(), Eigen::Vector3f{0, 0, 0});
        // 储存小像素的颜色信息
        for (int i = 0; i < frame_buf_2xSSAA.size(); i++) {
            frame_buf_2xSSAA[i].resize(4);
            std::fill(frame_buf_2xSSAA[i].begin(), frame_buf_2xSSAA[i].end(), Eigen::Vector3f{0, 0, 0});
        }
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth) {
        std::fill(depth_buf.begin(), depth_buf.end(), std::numeric_limits<float>::infinity());
        // 储存小像素的深度信息
        for (int i = 0; i < depth_buf_2xSSAA.size(); i++) {
            depth_buf_2xSSAA[i].resize(4);
            std::fill(depth_buf_2xSSAA[i].begin(), depth_buf_2xSSAA[i].end(), std::numeric_limits<float>::infinity());
        }
    }
}

rst::rasterizer::rasterizer(int h, int w) : height(h), width(w) {
    // frame_buf, depth_buf都是二维的，此处将他们都映射为一维，或许是为了方便openCV作图

    // 帧缓存，保存图像中每个像素显示什么颜色
    frame_buf.resize(w * h);

    // Z-buffer深度缓存，保存图像中每个像素在三维空间中的深度信息，后续处理遮挡关系
    // Depth-Buffer（深度缓存）有两种：Z-Buffer 和 W-Buffer
    // https://developer.aliyun.com/article/49272
    // https://www.csie.ntu.edu.tw/~r89004/hive/hsr/page_1.html
    depth_buf.resize(w * h);

    // SSAA超采样，每个像素划分成四个独立的小像素
    frame_buf_2xSSAA.resize(w * h);
    depth_buf_2xSSAA.resize(w * h);
}

int rst::rasterizer::get_index(int x, int y) {
    return (height - 1 - y) * width + x;
}

// 设置图像中某一个像素的颜色信息，光栅器的最后一步
void rst::rasterizer::set_pixel(const Eigen::Vector3f &point, const Eigen::Vector3f &color) {
    //old index: auto ind = point.y() + point.x() * width;
    auto ind = (height - 1 - point.y()) * width + point.x();
    frame_buf[ind] = color;
}

// clang-format on