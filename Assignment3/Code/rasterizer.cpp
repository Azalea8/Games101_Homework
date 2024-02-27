//
// Created by goksu on 4/6/19.
//

#include <algorithm>
#include "rasterizer.hpp"
#include <opencv2/opencv.hpp>
#include <cmath>

const bool SSAA = false;

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

rst::col_buf_id rst::rasterizer::load_normals(const std::vector<Eigen::Vector3f> &normals) {
    auto id = get_next_id();
    nor_buf.emplace(id, normals);

    normal_id = id;

    return {id};
}


// Bresenham's line drawing algorithm
void rst::rasterizer::draw_line(Eigen::Vector3f begin, Eigen::Vector3f end) {
    auto x1 = begin.x();
    auto y1 = begin.y();
    auto x2 = end.x();
    auto y2 = end.y();

    Eigen::Vector3f line_color = {255, 255, 255};

    int x, y, dx, dy, dx1, dy1, px, py, xe, ye, i;

    dx = x2 - x1;
    dy = y2 - y1;
    dx1 = fabs(dx);
    dy1 = fabs(dy);
    px = 2 * dy1 - dx1;
    py = 2 * dx1 - dy1;

    if (dy1 <= dx1) {
        if (dx >= 0) {
            x = x1;
            y = y1;
            xe = x2;
        } else {
            x = x2;
            y = y2;
            xe = x1;
        }
        Eigen::Vector2i point = Eigen::Vector2i(x, y);
        set_pixel(point, line_color);
        for (i = 0; x < xe; i++) {
            x = x + 1;
            if (px < 0) {
                px = px + 2 * dy1;
            } else {
                if ((dx < 0 && dy < 0) || (dx > 0 && dy > 0)) {
                    y = y + 1;
                } else {
                    y = y - 1;
                }
                px = px + 2 * (dy1 - dx1);
            }
//            delay(0);
            Eigen::Vector2i point = Eigen::Vector2i(x, y);
            set_pixel(point, line_color);
        }
    } else {
        if (dy >= 0) {
            x = x1;
            y = y1;
            ye = y2;
        } else {
            x = x2;
            y = y2;
            ye = y1;
        }
        Eigen::Vector2i point = Eigen::Vector2i(x, y);
        set_pixel(point, line_color);
        for (i = 0; y < ye; i++) {
            y = y + 1;
            if (py <= 0) {
                py = py + 2 * dx1;
            } else {
                if ((dx < 0 && dy < 0) || (dx > 0 && dy > 0)) {
                    x = x + 1;
                } else {
                    x = x - 1;
                }
                py = py + 2 * (dx1 - dy1);
            }
//            delay(0);
            Eigen::Vector2i point = Eigen::Vector2i(x, y);
            set_pixel(point, line_color);
        }
    }
}

auto to_vec4(const Eigen::Vector3f &v3, float w = 1.0f) {
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}

static bool insideTriangle(int x, int y, const Vector4f *_v) {
    Vector3f v[3];
    for (int i = 0; i < 3; i++)
        v[i] = {_v[i].x(), _v[i].y(), 1.0};
    Vector3f f0, f1, f2;
    f0 = v[1].cross(v[0]);
    f1 = v[2].cross(v[1]);
    f2 = v[0].cross(v[2]);
    Vector3f p(x, y, 1.);
    if ((p.dot(f0) * f0.dot(v[2]) > 0) && (p.dot(f1) * f1.dot(v[0]) > 0) && (p.dot(f2) * f2.dot(v[1]) > 0))
        return true;
    return false;
}

static std::tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector4f *v) {
    float c1 = (x * (v[1].y() - v[2].y()) + (v[2].x() - v[1].x()) * y + v[1].x() * v[2].y() - v[2].x() * v[1].y()) /
               (v[0].x() * (v[1].y() - v[2].y()) + (v[2].x() - v[1].x()) * v[0].y() + v[1].x() * v[2].y() -
                v[2].x() * v[1].y());
    float c2 = (x * (v[2].y() - v[0].y()) + (v[0].x() - v[2].x()) * y + v[2].x() * v[0].y() - v[0].x() * v[2].y()) /
               (v[1].x() * (v[2].y() - v[0].y()) + (v[0].x() - v[2].x()) * v[1].y() + v[2].x() * v[0].y() -
                v[0].x() * v[2].y());
    float c3 = (x * (v[0].y() - v[1].y()) + (v[1].x() - v[0].x()) * y + v[0].x() * v[1].y() - v[1].x() * v[0].y()) /
               (v[2].x() * (v[0].y() - v[1].y()) + (v[1].x() - v[0].x()) * v[2].y() + v[0].x() * v[1].y() -
                v[1].x() * v[0].y());
    return {c1, c2, c3};
}

void rst::rasterizer::draw(std::vector<Triangle *> &TriangleList) {

    float f1 = (50 - 0.1) / 2.0;
    float f2 = (-50 + -0.1) / 2.0;

    Eigen::Matrix4f mvp = projection * view * model;
    for (const auto &t: TriangleList) {
        Triangle newtri = *t;

        // 取得三维空间下的坐标 (还未投影到平面)
        std::array<Eigen::Vector4f, 3> mm{
                (view * model * t->v[0]),
                (view * model * t->v[1]),
                (view * model * t->v[2])
        };

        std::array<Eigen::Vector3f, 3> viewspace_pos;

        std::transform(mm.begin(), mm.end(), viewspace_pos.begin(), [](auto &v) {
            return v.template head<3>();
        });

        Eigen::Vector4f v[] = {
                mvp * t->v[0],
                mvp * t->v[1],
                mvp * t->v[2]
        };

        //Homogeneous division
        for (auto &vec: v) {
            vec.x() /= vec.w();
            vec.y() /= vec.w();
            vec.z() /= vec.w();
        }

        Eigen::Matrix4f inv_trans = (view * model).inverse().transpose();
        Eigen::Vector4f n[] = {
                inv_trans * to_vec4(t->normal[0], 0.0f),
                inv_trans * to_vec4(t->normal[1], 0.0f),
                inv_trans * to_vec4(t->normal[2], 0.0f)
        };

        //Viewport transformation
        for (auto &vert: v) {
            vert.x() = 0.5 * width * (vert.x() + 1.0);
            vert.y() = 0.5 * height * (vert.y() + 1.0);
            vert.z() = vert.z() * f1 + f2;
        }

        for (int i = 0; i < 3; ++i) {
            //screen space coordinates
            newtri.setVertex(i, v[i]);
        }

        for (int i = 0; i < 3; ++i) {
            //view space normal
            newtri.setNormal(i, n[i].head<3>());
        }

        newtri.setColor(0, 148,121.0,92.0);
        newtri.setColor(1, 148,121.0,92.0);
        newtri.setColor(2, 148,121.0,92.0);

        // setColor会将颜色归一化
//        newtri.setColor(0, 255, 0.0, 0.0);
//        newtri.setColor(1, 255, 0.0, 0.0);
//        newtri.setColor(2, 255, 0.0, 0.0);

        // Also pass view space vertice position
        rasterize_triangle(newtri, viewspace_pos);
    }

    if (SSAA) {
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                Eigen::Vector2i p = {(float) x, (float) y};
                Eigen::Vector3f color(0, 0, 0);
                for (int i = 0; i < 4; i++)
                    color += frame_buf_2xSSAA[get_index(x, y)][i];
                color /= 4;
                set_pixel(p, color);
            }
        }
    }
}

static Eigen::Vector3f
interpolate(float alpha, float beta, float gamma, const Eigen::Vector3f &vert1, const Eigen::Vector3f &vert2,
            const Eigen::Vector3f &vert3, float weight) {
    return (alpha * vert1 + beta * vert2 + gamma * vert3) / weight;
}

static Eigen::Vector2f
interpolate(float alpha, float beta, float gamma, const Eigen::Vector2f &vert1, const Eigen::Vector2f &vert2,
            const Eigen::Vector2f &vert3, float weight) {
    auto u = (alpha * vert1[0] + beta * vert2[0] + gamma * vert3[0]);
    auto v = (alpha * vert1[1] + beta * vert2[1] + gamma * vert3[1]);

    u /= weight;
    v /= weight;

    return Eigen::Vector2f(u, v);
}

//Screen space rasterization
void rst::rasterizer::rasterize_triangle(const Triangle &t, const std::array<Eigen::Vector3f, 3> &view_pos) {
    // TODO: From your HW3, get the triangle rasterization code.
    // TODO: Inside your rasterization loop:
    //    * v[i].w() is the vertex view space depth value z.
    //    * Z is interpolated view space depth for the current pixel
    //    * zp is depth between zNear and zFar, used for z-buffer

    // float Z = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
    // float zp = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
    // zp *= Z;

    // TODO: Interpolate the attributes:
    // auto interpolated_color
    // auto interpolated_normal
    // auto interpolated_texcoords
    // auto interpolated_shadingcoords

    // Use: fragment_shader_payload payload( interpolated_color, interpolated_normal.normalized(), interpolated_texcoords, texture ? &*texture : nullptr);
    // Use: payload.view_pos = interpolated_shadingcoords;
    // Use: Instead of passing the triangle's color directly to the frame buffer, pass the color to the shaders first to get the final color;
    // Use: auto pixel_color = fragment_shader(payload);
    auto v = t.toVector4();
    int min_x = INT_MAX;
    int max_x = INT_MIN;
    int min_y = INT_MAX;
    int max_y = INT_MIN;
    for (auto point: v) {
        if (point[0] < min_x) min_x = point[0];
        if (point[0] > max_x) max_x = ceil(point[0]);
        if (point[1] < min_y) min_y = point[1];
        if (point[1] > max_y) max_y = ceil(point[1]);
    }
    for (int x = min_x; x <= max_x; x++) {
        for (int y = min_y; y <= max_y; y++) {
            // SSAA抗锯齿后没有明显效果，运行时间从原来的 8秒变为 26秒
            if (SSAA) {
                int index = 0;
                // 划分四个小的像素
                for (float i = 0.25; i < 1.0; i += 0.5) {
                    for (float j = 0.25; j < 1.0; j += 0.5) {
                        if (insideTriangle(x + i, y + j, t.v)) {
                            //得到这个点的重心坐标
                            auto abg = computeBarycentric2D((float) x + i, (float) y + j, t.v);
                            float alpha = std::get<0>(abg);
                            float beta = std::get<1>(abg);
                            float gamma = std::get<2>(abg);
                            //z-buffer插值
                            float w_reciprocal = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w()); //归一化系数
                            float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() +
                                                   gamma * v[2].z() / v[2].w();
                            z_interpolated *= w_reciprocal;

                            if (abs(z_interpolated) < depth_buf_2xSSAA[get_index(x, y)][index]) {
                                Eigen::Vector2i p = {(float) x, (float) y};
                                // 颜色插值
                                auto interpolated_color = interpolate(alpha, beta, gamma, t.color[0], t.color[1],
                                                                      t.color[2], 1);
                                // 法向量插值
                                auto interpolated_normal = interpolate(alpha, beta, gamma, t.normal[0], t.normal[1],
                                                                       t.normal[2], 1);
                                // 纹理颜色插值
                                auto interpolated_texcoords = interpolate(alpha, beta, gamma, t.tex_coords[0],
                                                                          t.tex_coords[1], t.tex_coords[2], 1);
                                // 内部点位置插值
                                auto interpolated_shadingcoords = interpolate(alpha, beta, gamma, view_pos[0],
                                                                              view_pos[1], view_pos[2], 1);

                                fragment_shader_payload payload(interpolated_color, interpolated_normal.normalized(),
                                                                interpolated_texcoords, texture ? &*texture : nullptr);

                                payload.view_pos = interpolated_shadingcoords;
                                auto pixel_color = fragment_shader(payload);
                                // set_pixel(p, pixel_color); //设置颜色
                                frame_buf_2xSSAA[get_index(x, y)][index] = pixel_color;
                                depth_buf_2xSSAA[get_index(x, y)][index] = abs(z_interpolated);//更新z值
                            }
                        }
                        index++;
                    }
                }
            }
            else {
                //以像素中心点作为采样点
                if (insideTriangle((float) x + 0.5, (float) y + 0.5, t.v)) {

                    //得到这个点的重心坐标
                    auto abg = computeBarycentric2D((float) x + 0.5, (float) y + 0.5, t.v);
                    float alpha = std::get<0>(abg);
                    float beta = std::get<1>(abg);
                    float gamma = std::get<2>(abg);

                    //z-buffer插值
                    float w_reciprocal = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w()); //归一化系数
                    float z_interpolated =
                            alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                    z_interpolated *= w_reciprocal;

                    if (abs(z_interpolated) < depth_buf[get_index(x, y)]) {
                        Eigen::Vector2i p = {(float) x, (float) y};
                        // 颜色插值
                        auto interpolated_color = interpolate(alpha, beta, gamma, t.color[0], t.color[1],
                                                              t.color[2],1);
                        // 法向量插值
                        auto interpolated_normal = interpolate(alpha, beta, gamma, t.normal[0], t.normal[1],
                                                               t.normal[2],1);
                        // 纹理插值
                        auto interpolated_texcoords = interpolate(alpha, beta, gamma, t.tex_coords[0], t.tex_coords[1],
                                                                  t.tex_coords[2], 1);
                        // 内部点 3维空间插值
                        auto interpolated_shadingcoords = interpolate(alpha, beta, gamma, view_pos[0], view_pos[1],
                                                                      view_pos[2], 1);

                        fragment_shader_payload payload(interpolated_color, interpolated_normal.normalized(),
                                                        interpolated_texcoords, texture ? &*texture : nullptr);

                        payload.view_pos = interpolated_shadingcoords;
                        auto pixel_color = fragment_shader(payload);
                        set_pixel(p, pixel_color); //设置颜色
                        depth_buf[get_index(x, y)] = abs(z_interpolated);//更新z值
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

rst::rasterizer::rasterizer(int h, int w) : width(w), height(h) {
    frame_buf.resize(w * h);
    depth_buf.resize(w * h);

    // SSAA超采样，每个像素划分成四个独立的小像素
    frame_buf_2xSSAA.resize(w * h);
    depth_buf_2xSSAA.resize(w * h);

    texture = std::nullopt;
}

int rst::rasterizer::get_index(int x, int y) {
    return (height - 1 - y) * width + x;
}

void rst::rasterizer::set_pixel(const Vector2i &point, const Eigen::Vector3f &color) {
    //old index: auto ind = point.y() + point.x() * width;
    int ind = (height - 1 - point.y()) * width + point.x();
    frame_buf[ind] = color;
}

void rst::rasterizer::set_vertex_shader(std::function<Eigen::Vector3f(vertex_shader_payload)> vert_shader) {
    vertex_shader = vert_shader;
}

void rst::rasterizer::set_fragment_shader(std::function<Eigen::Vector3f(fragment_shader_payload)> frag_shader) {
    fragment_shader = frag_shader;
}

