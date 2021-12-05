# CMakeList.txt : CMake project for grid_test, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)


target_sources_local(grid_test 
	PRIVATE 
	   "ndgrid.h"
	   "grid3d.h"
	   "grid3d.cpp")