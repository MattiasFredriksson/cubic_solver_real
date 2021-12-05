# CMakeList.txt : CMake project for DTW, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Static library output: DTW
add_library(${PROJECT} STATIC "")

# Preprocessor options
target_compile_options(${PROJECT} PRIVATE -openmp)
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
target_compile_definitions(${PROJECT} PRIVATE CONTEXT_DEBUG_OUT)

target_compile_options(${PROJECT} PRIVATE -Zi)
else() # Release
target_compile_options(${PROJECT} PRIVATE -O2)
endif()

# Preprocessor defines
#target_compile_definitions(${PROJECT} PRIVATE "EIGEN_DEFAULT_TO_ROW_MAJOR")
 
# Header libs
#target_include_directories(${PROJECT} PRIVATE "../../../external/eigen")

# Subdirectories
add_subdirectory("include")
add_subdirectory("src")

# Inlude directories
target_include_directories(${PROJECT} PRIVATE "include")
