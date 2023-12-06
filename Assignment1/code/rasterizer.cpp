//
// Created by goksu on 4/6/19.
//

#include <algorithm>
#include "rasterizer.hpp"
#include <opencv2/opencv.hpp>
#include <math.h>
#include <stdexcept>


rst::pos_buf_id rst::rasterizer::load_positions(const std::vector<Eigen::Vector3f> &positions)
{
    auto id = get_next_id();
    pos_buf.emplace(id, positions); // 存入顶点缓存中

    return {id};
}

rst::ind_buf_id rst::rasterizer::load_indices(const std::vector<Eigen::Vector3i> &indices)
{
    auto id = get_next_id();
    ind_buf.emplace(id, indices);

    return {id};
}

// Bresenham's line drawing algorithm
// 离散划线算法，没有涉及到光栅化
// Code taken from a stack overflow answer: https://stackoverflow.com/a/16405254
void rst::rasterizer::draw_line(Eigen::Vector3f begin, Eigen::Vector3f end)
{
    auto x1 = begin.x();
    auto y1 = begin.y();
    auto x2 = end.x();
    auto y2 = end.y();

    Eigen::Vector3f line_color = {255, 255, 255};

    int x,y,dx,dy,dx1,dy1,px,py,xe,ye,i;

    dx=x2-x1;
    dy=y2-y1;
    dx1=fabs(dx);
    dy1=fabs(dy);
    px=2*dy1-dx1;
    py=2*dx1-dy1;

    if(dy1<=dx1)
    {
        if(dx>=0)
        {
            x=x1;
            y=y1;
            xe=x2;
        }
        else
        {
            x=x2;
            y=y2;
            xe=x1;
        }
        Eigen::Vector3f point = Eigen::Vector3f(x, y, 1.0f);
        set_pixel(point,line_color);
        for(i=0;x<xe;i++)
        {
            x=x+1;
            if(px<0)
            {
                px=px+2*dy1;
            }
            else
            {
                if((dx<0 && dy<0) || (dx>0 && dy>0))
                {
                    y=y+1;
                }
                else
                {
                    y=y-1;
                }
                px=px+2*(dy1-dx1);
            }
//            delay(0);
            Eigen::Vector3f point = Eigen::Vector3f(x, y, 1.0f);
            set_pixel(point,line_color);
        }
    }
    else
    {
        if(dy>=0)
        {
            x=x1;
            y=y1;
            ye=y2;
        }
        else
        {
            x=x2;
            y=y2;
            ye=y1;
        }
        Eigen::Vector3f point = Eigen::Vector3f(x, y, 1.0f);
        set_pixel(point,line_color);
        for(i=0;y<ye;i++)
        {
            y=y+1;
            if(py<=0)
            {
                py=py+2*dx1;
            }
            else
            {
                if((dx<0 && dy<0) || (dx>0 && dy>0))
                {
                    x=x+1;
                }
                else
                {
                    x=x-1;
                }
                py=py+2*(dx1-dy1);
            }
//            delay(0);
            Eigen::Vector3f point = Eigen::Vector3f(x, y, 1.0f);
            set_pixel(point,line_color);
        }
    }
}

// 齐次坐标
auto to_vec4(const Eigen::Vector3f& v3, float w = 1.0f)
{
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}

void rst::rasterizer::draw(rst::pos_buf_id pos_buffer, rst::ind_buf_id ind_buffer, rst::Primitive type)
{
    if (type != rst::Primitive::Triangle) // 只接受三角形类型
    {
        throw std::runtime_error("Drawing primitives other than triangle is not implemented yet!");
    }

    // 从顶点缓存中读取当前帧所需要绘制的三角形信息

    // pos_id为唯一标识, pos_buffer为结构体封装 pos_id, pos_buf作为 Map，对应光栅器缓存的三角形数据
    auto& buf = pos_buf[pos_buffer.pos_id];
    // 与顶点坐标类似，用 ind来分组构造有序的三角形
    auto& ind = ind_buf[ind_buffer.ind_id];

    float f1 = (100 - 0.1) / 2.0;
    float f2 = (100 + 0.1) / 2.0;

    Eigen::Matrix4f mvp = projection * view * model; // 做 MVP变换，矩阵运算具有结合律，所以可以先算

    for (auto& i : ind) // 三角形原先分过组，每组一次变换，一一对应
    {
        Triangle t; // 实例化一个三角形类

        // v是个数据类型为4维向量(点的齐次坐标)的数组, 共有三个元素代表三角形的顶点。 v数组对应MVP变换后的三角形
        Eigen::Vector4f v[] = {
                mvp * to_vec4(buf[i[0]], 1.0f),
                mvp * to_vec4(buf[i[1]], 1.0f),
                mvp * to_vec4(buf[i[2]], 1.0f),
        };

        // 对齐次坐标归一化
        for (auto& vec : v) {
            vec /= vec.w();
        }

        // 视口变换 ([-1, 1]^2变换到显示器)
        for (auto & vert : v)
        {
            vert.x() = 0.5 * width * (vert.x() + 1.0);
            vert.y() = 0.5 * height * (vert.y() + 1.0);
            vert.z() = vert.z() * f1 + f2; // Z为深度信息, 这里的变换我猜测是防止Z值太小, 对后面的遮挡判断有影响
        }

        // 将变换后的三角形交给三角形类
        for (int j = 0; j < 3; ++j)
        {
            t.setVertex(j, v[j].head<3>());
            // t.setVertex(i, v[i].head<3>());
            // t.setVertex(i, v[i].head<3>());
            // 不知道为什么要三次同样的的操作
        }

        // 设定当前三角形的顶点颜色
        t.setColor(0, 120.0,  0.0,  0.0);
        t.setColor(1, 0.0  ,120.0,  0.0);
        t.setColor(2, 0.0  ,  0.0,120.0);

        // 绘制当前的三角形的三边
        rasterize_wireframe(t);
    }
}

void rst::rasterizer::rasterize_wireframe(const Triangle& t)
{
    // 画线， 从三角形的顶点到顶点
    draw_line(t.c(), t.a());
    draw_line(t.c(), t.b());
    draw_line(t.b(), t.a());
}

void rst::rasterizer::set_model(const Eigen::Matrix4f& m)
{
    model = m;
}

void rst::rasterizer::set_view(const Eigen::Matrix4f& v)
{
    view = v;
}

void rst::rasterizer::set_projection(const Eigen::Matrix4f& p)
{
    projection = p;
}

void rst::rasterizer::clear(rst::Buffers buff)
{
    if ((buff & rst::Buffers::Color) == rst::Buffers::Color)
    {
        std::fill(frame_buf.begin(), frame_buf.end(), Eigen::Vector3f{0, 0, 0});
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth)
    {
        std::fill(depth_buf.begin(), depth_buf.end(), std::numeric_limits<float>::infinity());
    }
}

rst::rasterizer::rasterizer(int w, int h) : width(w), height(h)
{
    frame_buf.resize(w * h); //设置帧缓存大小
    depth_buf.resize(w * h); //设置深度缓存
}

int rst::rasterizer::get_index(int x, int y)
{
    // opencv的坐标是传统的二维数组，这里我们的像素是高中传统的坐标系
    return (height-y)*width + x;
}

void rst::rasterizer::set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color)
{
    //old index: auto ind = point.y() + point.x() * width;
    if (point.x() < 0 || point.x() >= width || point.y() < 0 || point.y() >= height){
        return;
    }
    auto ind = (height-point.y())*width + point.x();
    frame_buf[ind] = color;
}

