# tmc

**Construction of Topologically Correct and Manifold Iso-surfaces**
Author: Roberto Grosso

Different algorithms for computing topologically correct and manifold iso-surfaces from volume data (hexahedral meshes or voxel grids) were implemented. The output surfaces are represented using an indexed face set (shared vertex) and a halfedge data structure. The acronym **tmc** is intended to mean *topologically correct and manifold iso-surface extraction by using a Marching Cubes like algorithm*.

## Folder matlab

In this folder there are some few MATLAB scripts which computes the intersection of a level set with one or two unit cells for different configurations.

### Script t_mc_figs.m

This script generates the figures presented in the paper [1]. Just open the MATLAB script and press run. It will generate a set of  MATLAB `*.fig` figures. 

### Scripts t_mc_gui.m and t_mc.m

These scripts compute the intersection of a level set with a unit cell for different configurations as described in the paper [1]. They contain a user interface to facilitate the usage and visualize different cases for different configurations. In order to run the script use MATLAB and run the script **tmc_gui.m**.

### Scripts singular_gui.m and t_singular.m

These scripts implement the *asymptotic decider* presented in [2]. Some special cases handled by the asymptotic decider are demonstrated in this script. They can be executed by running the user interface in **singular_gui.m**

## Folder mc

In this folder you find the main function in tmc.cpp and the classes in MarchingCubes.h and MarchingCubes.cpp which compute the topologically correct iso-surface from volume data as presented in [1]. Furthermore, a small data set used in the paper. The file format used for the volume data can be easily read from the method which overload the `function operator`. A shell script is provided to compile and execute the code.

Execution:  

1) within the folder containing the sources just run the bash-script `runtmc.sh`.

2) Output are the files `skull128.obj` and `skull128.off`.

## Folder cuda_tmc

This folder contains a CUDA implementation of the topologically correct iso-surface algorithm presented in [1]. This implementation optionally generates a halfedge or a shared vertex data structure for the triangle mesh. We do not include a `main` function to run the code. The class `MarchingCubed` implements two different methods to generate iso-surfaces from volume data: `mc_halfedge()` generates a halfedge data structure and `mc_sharedvertex()` generates a indexed face set (shared vertex) data structure.

## Folder dmc

This folder contains a CUDA implementation of the *Dual Marching Cubes* method presented in [3]. The main function demonstrates how to use the code. It can read volume data or alternatively generate a synthetic data set for demonstration purposes. The file format for volume data can be easily seen in the method `readDataFromFile()` implemented in the class `DualMarchingCubes`.

[1]: Roberto Grosso: **Construction of Topologically Correct and Manifold Isosurfaces**. *Computer Graphics Forum 35(5):187-196 · August 2016*

[2]: Roberto Grosso: **An asymptotic decider for robust and topologically correct triangulation of isosurfaces: topologically correct isosurfaces**. *CGI '17 Proceedings of the Computer Graphics International Conference*. Japan, June 2017

[3]: Roberto Grosso, Daniel Zint: **PARALLEL RECONSTRUCTION OF QUAD ONLY MESHES FROM VOLUME DATA**. *Proceedings of the 28th International Meshing Roundtable*. Buffalo, NY, USA, October 2019.
