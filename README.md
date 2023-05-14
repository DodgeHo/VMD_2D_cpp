# VMD_2D_cpp
# C++ implementation of 2D Variational Mode Decomposition using Eigen3
Written by: Dodge(Lang HE) asdsay@gmail.com
Updated date: 2023-05-14

VMD, aka Variational Mode Decomposition, is a signal processing tool that decompse the input signal into different band-limited IMFs. 
Similar to the Project [**VMD_CPP**](https://github.com/DodgeHo/VMD_cpp), Project **VMD_2D_CPP**  is an imitation of [that in MATLAB](https://www.mathworks.com/matlabcentral/fileexchange/45918-two-dimensional-variational-mode-decomposition). Spectrum-based decomposition of a 2D input signal into k band-separated modes. 

If you are looking for document to describe Variational mode decomposition, please turn to the original paper [Variational Mode Decomposition](ftp://ftp.math.ucla.edu/pub/camreport/cam14-16.pdf). You can also find the MATLAB codes here.

In this project, I used Eigen3 to refactor VMD in C++, so that we can use it without MATLAB. Also here, I input an image by CImg to test, but you don't need to install CImg, only Eigen. This sample code was written in MSBuild. You can both run in Visual Studio 2022 or MSVC or CMAKE/GCC, you can use either the sln project file, or the CMakeList.txt, they both work.
Detail input and output please check out function **VMD_2D** in file VMD_2D_Utils.cpp.

PS: This VMD_2D runs too slow. I will use OpenCV to re-write it soon.
