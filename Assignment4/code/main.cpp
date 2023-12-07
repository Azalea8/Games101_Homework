#include <chrono>
#include <iostream>
#include <opencv2/opencv.hpp>

std::vector<cv::Point2f> control_points;

// 贝塞尔曲线的控制点数量，默认为 4，即 3阶贝塞尔曲线
#define NUM 4

void mouse_handler(int event, int x, int y, int flags, void *userdata) 
{
    if (event == cv::EVENT_LBUTTONDOWN && control_points.size() < NUM)
    {
        std::cout << "Left button of the mouse is clicked - position (" << x << ", " << y << ")" << '\n';
        control_points.emplace_back(x, y);
    }     
}

// 显式调用，没有递归
void naive_bezier(const std::vector<cv::Point2f> &points, cv::Mat &window) 
{
    auto &p_0 = points[0];
    auto &p_1 = points[1];
    auto &p_2 = points[2];
    auto &p_3 = points[3];

    for (double t = 0.0; t <= 1.0; t += 0.001) 
    {
        auto point = std::pow(1 - t, 3) * p_0 + 3 * t * std::pow(1 - t, 2) * p_1 +
                 3 * std::pow(t, 2) * (1 - t) * p_2 + std::pow(t, 3) * p_3;

        // 曲线颜色设置为红色，OpenCV为 BGR
        window.at<cv::Vec3b>(point.y, point.x)[2] = 255;
    }
}

// 贝塞尔曲线的基本公式
cv::Point2f lerp_v2f(const cv::Point2f& a, const cv::Point2f& b, float t)
{
    return (1 - t) * a + t * b;
}

// 递归
cv::Point2f recursive_bezier(const std::vector<cv::Point2f>& control_points, float t)
{
    // TODO: Implement de Casteljau's algorithm
    // 递归出口，此时为贝塞尔曲线上的点
    if (control_points.size() == 1)
    {
        return control_points[0];
    }

    // 新建一个数组，储存新的控制点，数量每次会减少 1，减少到 1时即为出口
    std::vector<cv::Point2f> lerp_points;
    // 遍历原控制点数组，将通过 lerp_v2f（最基本的公式）新形成的控制点依次填入新的数组
    for (size_t i = 1;i < control_points.size();i++)
    {
        lerp_points.push_back(lerp_v2f(control_points[i - 1], control_points[i], t));
    }
    // 新数组继续递归
    return recursive_bezier(lerp_points, t);
}

// 遍历参数 t
void bezier(const std::vector<cv::Point2f> &control_points, cv::Mat &window) 
{
    // TODO: Iterate through all t = 0 to t = 1 with small steps, and call de Casteljau's 
    // recursive Bezier algorithm.
    for (double t = 0.0;t <= 1.0; t += 0.001)
    {
        // 对于每一个 t，求它对应的贝塞尔曲线的点，这里是通过最基本的公式 (1 - t) * a + t * b，递归求得，
        // 实际上也可以直接用 n阶展开式，直接求得，就像 naive_bezier函数那样
        auto point = recursive_bezier(control_points, t);

        // 坐标系转换，吐槽一句，opencv鼠标点击 和 opencv储存图像的坐标系不一样
        window.at<cv::Vec3b>(point.y, point.x)[2] = 255;
    }
}

// 计算阶乘
int factorial(int n)
{
    int fc=1;
    for(int i=1;i<=n;++i) fc *= i;
    return fc;
}

// 计算组合数
/*
    从 n个不同元素中取出 m ( m ≤ n )个元素的所有组合的个数，叫做 n个不同元素中取出 m个元素的组合数。用符号 c(m,n) 表示。
    组合数公式：c(m,n) = n! / (m! * (n-m)!)
    性质：c(n,m) = c(n,m-n)
    递推公式：c(n,m) = c(n-1,m-1) + c(n-1,m)
*/
int combo(int m,int n)
{
    int com = factorial(n) / (factorial(m) * factorial(n-m));
    return com;
}

void bezier_combo(const std::vector<cv::Point2f> &points, cv::Mat &window) {
    for (double t = 0.0; t <= 1.0; t += 0.001) {
        cv::Point2f point(0, 0);
        for (size_t i = 0;i < points.size();i++) {
            point += combo(i, NUM - 1) * std::pow((1 - t), NUM - 1 - i) * std::pow(t, i) * points[i];
        }
        // 曲线颜色设置为红色，OpenCV为 BGR
        window.at<cv::Vec3b>(point.y, point.x)[2] = 255;
    }
}

int main() 
{
    cv::Mat window = cv::Mat(700, 700, CV_8UC3, cv::Scalar(0));
    cv::cvtColor(window, window, cv::COLOR_BGR2RGB);
    cv::namedWindow("Bezier Curve", cv::WINDOW_AUTOSIZE);

    cv::setMouseCallback("Bezier Curve", mouse_handler, nullptr);

    int key = -1;
    while (key != 27) 
    {
        for (auto &point : control_points) 
        {
            cv::circle(window, point, 3, {255, 255, 255}, 3);
        }

        if (control_points.size() == NUM)
        {
            // naive_bezier(control_points, window);
            // bezier(control_points, window);
            bezier_combo(control_points, window);

            cv::imshow("Bezier Curve", window);
            cv::imwrite("my_bezier_curve.png", window);
            key = cv::waitKey(0);

            return 0;
        }

        cv::imshow("Bezier Curve", window);
        key = cv::waitKey(20);
    }

return 0;
}
