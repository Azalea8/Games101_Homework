# Games101_Homework

### 所有编译好的依赖都放在了 *Games101_environment* 文件夹下
* [课程地址](https://www.bilibili.com/video/BV1X7411F744)
* [框架解析](https://space.bilibili.com/523296472)
---

## OpenCV
由于我采用的 *MinGW* 作为编译器，需要 *Cmake* 编译 *OpenCV* 的源代码，由于种种原因，可能用 *Cmake* 编译会出现各种问题，万幸的是，有人已经上传了 *OpenCV-MinGW-Build* 各种版本

* 下载地址 [*OpenCV-MinGW-Build*](https://github.com/huihut/OpenCV-MinGW-Build)

该项目所用的版本号为 *4.5.2*

## Eigen3
* 下载地址 [*Eigen*](https://eigen.tuxfamily.org/index.php?title=Main_Page)

* 解压到任意目录，*eigen* 的根目录下新建一个 *build* 目录

* 打开 *cmake-gui*，*source* 目录设置成  *eigen* 的根目录，*build* 目录设置成刚刚新建的目录，然后点击 *configure*，出来的 *makefile* 的格式选择 *mingw*，等待 *config* 结束

* *install* 地址也可以改成其他路径，尽量不要放在 *C* 盘，避免日后找不到

* 点击 *generate*

* 管理员权限运行 *cmd*，进入 *build* 目录，运行 *mingw32-make*,之后运行 *mingw32-make install*

* 删除解压出来的这个 *eigen* 目录(因为上面安装的已经安装到了 *c* 盘 *Program File x86* 里了，有时候 *build* 的时候会找错路径)

#### 记得在 *CMakeLists.txt* 中修改路径

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

## FreeType
* 下载地址 [*freetype*](https://download.savannah.gnu.org/releases/freetype/)

* 解压到任意目录，*eigen* 的根目录下新建一个 *build* 目录

* 打开 *cmake-gui*，*source* 目录设置成  *eigen* 的根目录，*build* 目录设置成刚刚新建的目录，然后点击 *configure*，出来的 *makefile*的格式选择 *mingw*，等待 *config*结束

* 然后会有很多选项，把能勾的都勾上，*install* 地址也可以改成其他路径，尽量不要放在 *C*盘，避免日后找不到

* 点击 *generate*

* 管理员权限运行 *cmd*，进入 *build* 目录，运行 *mingw32-make*,之后运行 *mingw32-make install*

* 这时候就会在你设置好的目录下看到如下项目结构

        include
            freetype2/..
        lib
            cmake/..
            pkgconfig/..
            libfreetype.a

#### 作业八改动的地方相对较多，这里就不贴具体改动的部分，可以直接用我上传的 *Assignment8* 充当作业的框架

