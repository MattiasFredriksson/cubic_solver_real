# CMakeList.txt : CMake project for grid_test, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Add source to this project's executable.
target_sources_local(${PROJECT_CTEST} 
	PRIVATE 
	   "main.cpp")