# tmc

Paper:


Author
Roberto Grosso

Platform:
The code was tested on a Mac Pro, with OS X El Capitan Version 10.11.5 and MATLAB 2016a and 
clang++ Apple LLVM version 7.3.0 (clang-703.0.29)

Folder matlab
Within this folder run the following command in order to generate MATLAB fig containing the data used to generate the 
figures in the paper.

In this folder you find the MATLAB script t_mc_figs.m wich compute the intersection of the level set with a unit cell.


Folder mc
Here you find the C++ code. Run the bash-script runtmc.sh to compile and run the code. The output are the files 
"skull128.obj" and "skull128.off". These data sets were used to generate the figures in the paper.

In this folder you find the main function in tmc.cpp and the classes in MarchingCubes.h and MarchingCubes.cpp which compute
the topologically correct isosurface from an uniform grid. Furthermore, I provide a small data set used in the paper. Larger
data sets were omitted due to their file sizes and can be delivered at any time if required.
