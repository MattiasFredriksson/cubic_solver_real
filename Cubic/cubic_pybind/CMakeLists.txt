# CMakeList.txt : CMake project for DTW_python, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Python dll output: ${PROJECT}_python
add_subdirectory ("pybind11")
pybind11_add_module(${PROJECT_PYTHON} "")
# Project name definition
target_compile_definitions(${PROJECT_PYTHON} PRIVATE PROJECT_NAME_DEF=${PROJECT_PYTHON})

# Include & link static library DTW
#find_package(pybind11)
target_include_directories(${PROJECT_PYTHON} PRIVATE "../${PROJECT_SDIR}/include")
target_link_libraries(${PROJECT_PYTHON} PRIVATE ${PROJECT})

# Header libs
#target_include_directories(${PROJECT_PYTHON} PRIVATE "../../../external/eigen")

# Preprocessor defines
target_compile_definitions(${PROJECT_PYTHON} PRIVATE "")
# Debug
#target_compile_definitions(${PROJECT_PYTHON} PRIVATE DEBUG_OUT)

# Preprocessor options
#target_compile_options(${PROJECT_PYTHON} PRIVATE -openmp)
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
target_compile_definitions(${PROJECT} PRIVATE CONTEXT_DEBUG_OUT)

target_compile_options(${PROJECT} PRIVATE -Zi)
else() # Release
target_compile_options(${PROJECT} PRIVATE -O2)
endif()

# Compilation options
add_compile_options($<$<CXX_COMPILER_ID:MSVC>:/MP>)

# Sources
target_sources(${PROJECT_PYTHON} 
	PRIVATE 
	   "main.cpp")