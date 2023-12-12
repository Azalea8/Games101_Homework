# Games101_Homework

### 该项目的前3个作业需要 OpenCV 和 Eigen3 的第三方库

## OpenCV
由于我采用的 **gcc** 作为编译器，需要 **Cmake** 编译 **OpenCV** 的源代码，由于种种原因，可能用 **Cmake** 编译会出现各种问题，万幸的是，有人已经上传了 **OpenCV-MinGW-Build** 各种版本

* 下载地址 [OpenCV-MinGW-Build](https://github.com/huihut/OpenCV-MinGW-Build)

该项目所用的版本号为 **4.5.2**

记得在 **CMakeLists.txt** 中修改路径

    cmake_minimum_required(VERSION 3.10)
    project(Rasterizer)

    # 这里的路径记得改成你自己的路径
    set(OpenCV_DIR D:/OpenCV-MinGW-Build-OpenCV-4.5.2-x64) 
    find_package(OpenCV REQUIRED)

    set(CMAKE_CXX_STANDARD 17)

    # mingw32-make install命令执行后安装到了c盘Program File x86
    # 大概率是一样的路径
    include_directories("C:/Program Files (x86)/Eigen3/include/eigen3")

    add_executable(Rasterizer main.cpp rasterizer.hpp rasterizer.cpp Triangle.hpp Triangle.cpp)
    target_link_libraries(Rasterizer ${OpenCV_LIBRARIES})
    

## Eigen3

#### Eigen库编译
* 下载地址 [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)

* 解压到任意目录，eigen的根目录下新建一个build目录

* 打开cmake软件，source目录设置成eigen的根目录，build目录设置成刚刚新建的目录，然后点击configure，出来的makefile的格式选择mingw，等待config结束

* 点击generate

* 管理员权限运行cmd，进入build目录，运行mingw32-make,之后运行mingw32-make install

* 删除解压出来的这个eigen目录(因为上面安装的已经安装到了c盘Program File x86里了，有时候build的时候会找错路径)

