#pragma once

#include<array>
#include<cstdint>

template <typename FP>
struct Vec3 {
	FP x, y, z;
};

template <typename FP = double>
FP* grid3d(std::array<std::int64_t, 3> shape, std::array<FP, 3> size);

template <typename FP = double>
Vec3<FP>* grid3d_struct(std::array<std::int64_t, 3> shape, std::array<FP, 3> size);

template <typename FP>
Vec3<FP>* grid3d_struct_cpy(std::array<std::int64_t, 3> shape, std::array<FP, 3> size);

template <typename FP = double>
FP* grid3d_r(std::array<std::int64_t, 3> shape, std::array<FP, 3> size);

template <typename FP = double>
FP* grid3d_b(std::array<std::int64_t, 3> shape, std::array<FP, 3> size);

template <typename FP = double>
FP* grid3d_p1(std::array<std::int64_t, 3> shape, std::array<FP, 3> size);

template <typename FP = double>
FP* grid3d_p2(std::array<std::int64_t, 3> shape, std::array<FP, 3> size);


template <typename FP = double>
FP* grid3d_ins(std::array<std::int64_t, 3> shape, std::array<FP, 3> size);


template <typename FP>
FP* grid3d_memcpy(std::array<std::int64_t, 3> shape, std::array<FP, 3> size);