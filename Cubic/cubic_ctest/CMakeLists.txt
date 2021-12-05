# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)



###########
# Target(s)
###########
add_executable(${PROJECT_CTEST} "")

# Compiler options
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
else()
target_compile_options(${PROJECT_CTEST} PRIVATE -O2)
endif()
target_compile_options(${PROJECT_CTEST} PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/MP>)
#target_compile_options(${PROJECT_CTEST} PRIVATE -openmp)

target_include_directories(${PROJECT_CTEST} PRIVATE "../cubic_lib/include/")
target_link_libraries(${PROJECT_CTEST} PRIVATE ${PROJECT})

##########
# External
##########


##########

# Include project src files.
add_subdirectory ("src")
