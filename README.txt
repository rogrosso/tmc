# tmc

Paper: paper1074
Title: Construction of Topologically Correct and Manifold Isosurfaces
    
Author: Roberto Grosso

Platform:
The code was tested on a Mac Pro, with OS X El Capitan Version 10.11.5 and MATLAB 2016a and 
clang++ Apple LLVM version 7.3.0 (clang-703.0.29)

Folder matlab
In this folder you find the MATLAB script t_mc_figs.m wich compute the intersection of the level set with a unit cell.
Execution:
The MATLAB script should be run as follows:
1) Change to the folder where you have saved the MATLAB script t_mc_figs.m
2) run the command:
        /Applications/MATLAB_R2016a.app/bin/matlab -nodisplay -nosplash -nodesktop -r "run t_mc_figs.m ; exit"

Folder mc
In this folder you find the main function in tmc.cpp and the classes in MarchingCubes.h and MarchingCubes.cpp which compute
the topologically correct isosurface from an uniform grid. Furthermore, I provide a small data set used in the paper. Larger
data sets were omitted due to their file sizes and can be delivered at any time if required.

Execution:  
1) within this folder just run bash-script runtmc.sh to compile and run the code.
2) The output are the files "skull128.obj" and "skull128.off". These data sets were used to generate the figures in the paper.