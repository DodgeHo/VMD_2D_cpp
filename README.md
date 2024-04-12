# VMD_2D_cpp
# C++ implementation of 2D Variational Mode Decomposition using Eigen3
Written by: Dodge(Lang HE) asdsay@gmail.com
Updated date: 2023-11-15

VMD, aka Variational Mode Decomposition, is a signal processing tool that decompse the input signal into different band-limited IMFs. 
Similar to the Project [**VMD_CPP**](https://github.com/DodgeHo/VMD_cpp), Project **VMD_2D_CPP**  is an imitation of [that in MATLAB](https://www.mathworks.com/matlabcentral/fileexchange/45918-two-dimensional-variational-mode-decomposition). Spectrum-based decomposition of a 2D input signal into k band-separated modes. 

If you are looking for document to describe Variational mode decomposition, please turn to the original paper [Variational Mode Decomposition](ftp://ftp.math.ucla.edu/pub/camreport/cam14-16.pdf). You can also find the MATLAB codes here.

In this project, I used Eigen3 to refactor VMD in C++, so that we can use it without MATLAB. Also here, I input an image by CImg to test. This sample code was written in MSBuild. You can both run in Visual Studio 2022 or MSVC or CMAKE/GCC, you can use either the sln project file, or the CMakeList.txt, they both work.
Detail input and output please check out function **VMD_2D** in file VMD_2D_Utils.cpp.

PS: This VMD_2D runs too slow. I will use OpenCV to re-write it soon.

# 二维VMD（变分模态分解）的C++实现，使用了Eigen3

作者：Dodge asdsay@gmail.com 
更新日期：2023-11-15

VMD（变分模态分解）是一种信号处理算法，可以将输入信号分解为不同带限的内禀模态函数（IMFs）。类似于项目[**VMD_CPP**](https://github.com/DodgeHo/VMD_cpp)，本项目**VMD_2D_CPP**是参考于[其在MATLAB中的实现](https://www.mathworks.com/matlabcentral/fileexchange/45918-two-dimensional-variational-mode-decomposition)。基于频谱的二维输入信号分解为k个带分离模式。

如果需要描述变分模态分解的文档，可参阅原始论文[Variational Mode Decomposition](ftp://ftp.math.ucla.edu/pub/camreport/cam14-16.pdf)。此链接中有MATLAB代码。

在这个项目中，我用Eigen3实现了C++中的VMD，无须MATLAB。同时，为了输入一张图片进行测试，使用了另一个库CImg。这个示例代码是用MSBuild编写的。
本项目用MSBuild编写，也可在Visual Studio 2022、MSVC或CMAKE/GCC中运行，sln项目文件或CMakeList.txt都能用。详细的输入和输出请查看[VMD_2D_Utils.cpp](https://github.com/DodgeHo/VMD_2D_cpp/blob/main/VMD_2D_Utils.cpp)文件中的**VMD_2D_CPP**函数。

附：这个VMD_2D运行稍慢。我之后会尝试用OpenCV重新实现。
