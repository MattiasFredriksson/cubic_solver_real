#include"grid3d.h"
#include<assert.h>
#include "IndexRange.h"

template <typename FP>
FP* grid3d(std::array<std::int64_t, 3> shape, std::array<FP, 3> size) {
	/* Reference implementation. */
	using namespace std;
	/* Initialization. */
	int64_t nelem = 1;
	std::array<FP, 3> delta;
	std::array<FP, 3> shift;
	for (int64_t i = 0; i < 3; i++) {
		nelem *= shape[i];
		int64_t shn1 = shape[i] - 1;
		delta[i] = shn1 == 0 ? 0.0f : size[i] / shn1;
		shift[i] = shn1 / FP(2);
	}
	FP* data = new FP[nelem * 3];


	/* Fill arrays */
	FP* dptr = data;
	for (int64_t x = 0; x < shape[0]; x++) {
		FP xval = (x - shift[0]) * delta[0];
		for (int64_t y = 0; y < shape[1]; y++) {
			FP yval = (y - shift[1]) * delta[1];
			for (int64_t z = 0; z < shape[2]; z++) {
				*dptr++ = xval;
				*dptr++ = yval;
				*dptr++ = (z - shift[2]) * delta[2];
			}
		}
	}
	return data;
}
template float* grid3d(std::array<std::int64_t, 3> shape, std::array<float, 3> size);


template <typename FP>
Vec3<FP>* grid3d_struct(std::array<std::int64_t, 3> shape, std::array<FP, 3> size) {
	/* Struct, ~identical. */
	using namespace std;
	/* Initialization. */
	int64_t nelem = 1;
	std::array<FP, 3> delta;
	std::array<FP, 3> shift;
	for (int64_t i = 0; i < 3; i++) {
		nelem *= shape[i];
		int64_t shn1 = shape[i] - 1;
		delta[i] = shn1 == 0 ? 0.0f : size[i] / shn1;
		shift[i] = shn1 / FP(2);
	}
	Vec3<FP>* data = new Vec3<FP>[nelem];


	/* Fill arrays */
	//Vec3<FP>* dptr = data;
	for (int64_t x = 0; x < shape[0]; x++) {
		Vec3<FP>* dptr = data + (x * shape[1] * shape[2]);
		FP xval = (x - shift[0]) * delta[0];
		for (int64_t y = 0; y < shape[1]; y++) {
			FP yval = (y - shift[1]) * delta[1];
			for (int64_t z = 0; z < shape[2]; z++) {
				*dptr++ = Vec3<FP>{ xval, yval , (z - shift[2]) * delta[2] };
			}
		}
	}
	return data;
}
template Vec3<float>* grid3d_struct(std::array<std::int64_t, 3> shape, std::array<float, 3> size);


template <typename FP>
FP* grid3d_memcpy(std::array<std::int64_t, 3> shape, std::array<FP, 3> size) {
	/* Optimal single thread performance (reduced parallel gain). */
	using namespace std;
	/* Initialization. */
	int64_t nelem = 1;
	std::array<FP, 3> delta;
	std::array<FP, 3> shift;
	for (int64_t i = 0; i < 3; i++) {
		nelem *= shape[i];
		int64_t shn1 = shape[i] - 1;
		delta[i] = shn1 == 0 ? 0.0f : size[i] / shn1;
		shift[i] = shn1 / FP(2);
	}
	FP* data = new FP[nelem * 3];

	/* Fill first 2D slice */
	FP* dptr = data;
	FP zval = (0 - shift[2]) * delta[2];
	for (int64_t y = 0; y < shape[1]; y++) {
		for (int64_t x = 0; x < shape[2]; x++) {
			*dptr++ = (x - shift[0]) * delta[0];
			*dptr++ = (y - shift[1]) * delta[1];
			*dptr++ = zval;
		}
	}

	/* Fill remaining slices */
	const int64_t nfp_2d = shape[1] * shape[2] * 3;
	const int64_t cpy_size = nfp_2d * sizeof(FP);
	for (int64_t z = 1; z < shape[2]; z++) {
		FP* dptr = data + (z * nfp_2d);
		memcpy(reinterpret_cast<void*>(dptr), reinterpret_cast<void*>(data), cpy_size);
		FP zval = (z - shift[2]) * delta[2];
		for (; dptr < data + (z * (nfp_2d + 1)); dptr += 3) {
			*dptr = zval;
		}
	}
	return data;
}
template float* grid3d_memcpy(std::array<std::int64_t, 3> shape, std::array<float, 3> size);

template <typename FP>
FP* grid3d_p1(std::array<std::int64_t, 3> shape, std::array<FP, 3> size) {
	using namespace std;
	/* Initialization. */
	int64_t nelem = 1;
	std::array<FP, 3> delta;
	std::array<FP, 3> shift;
	for (int64_t i = 0; i < 3; i++) {
		nelem *= shape[i];
		int64_t shn1 = shape[i] - 1;
		delta[i] = shn1 == 0 ? 0.0f : size[i] / shn1;
		shift[i] = shn1 / FP(2);
	}
	FP* data = new FP[nelem * 3];


	/* Fill arrays */
#pragma omp parallel for
	for (int64_t x = 0; x < shape[0]; x++) {
		FP* dptr = data + (x * shape[1] * shape[2] * 3);
		FP xval = (x - shift[0]) * delta[0];
		for (int64_t y = 0; y < shape[1]; y++) {
			FP yval = (y - shift[1]) * delta[1];
			for (int64_t z = 0; z < shape[2]; z++) {
				*dptr++ = xval;
				*dptr++ = yval;
				*dptr++ = (z - shift[2]) * delta[2];
			}
		}
	}
	return data;
}
template float* grid3d_p1(std::array<std::int64_t, 3> shape, std::array<float, 3> size);

template <typename FP>
FP* grid3d_p2(std::array<std::int64_t, 3> shape, std::array<FP, 3> size) {
	/* Parallel second loop. */
	using namespace std;
	/* Initialization. */
	int64_t nelem = 1;
	std::array<FP, 3> delta;
	std::array<FP, 3> shift;
	for (int64_t i = 0; i < 3; i++) {
		nelem *= shape[i];
		int64_t shn1 = shape[i] - 1;
		delta[i] = shn1 == 0 ? 0.0f : size[i] / shn1;
		shift[i] = shn1 / FP(2);
	}
	FP* data = new FP[nelem * 3];


	/* Fill arrays */
	for (int64_t x = 0; x < shape[0]; x++) {
#pragma omp parallel for
		for (int64_t y = 0; y < shape[1]; y++) {
			FP* dptr = data + (x * shape[1] * shape[2] + y * shape[2]);
			for (int64_t z = 0; z < shape[2]; z++) {
				*dptr++ = (x - shift[0]) * delta[0];
				*dptr++ = (y - shift[1]) * delta[1];
				*dptr++ = (z - shift[2]) * delta[2];
			}
		}
	}
	return data;
}
template float* grid3d_p2(std::array<std::int64_t, 3> shape, std::array<float, 3> size);


template <typename FP>
Vec3<FP>* grid3d_struct_cpy(std::array<std::int64_t, 3> shape, std::array<FP, 3> size) {
	/* Struct + memcpy, ~perf. as struct. */
	using namespace std;
	/* Initialization. */
	int64_t nelem = 1;
	std::array<FP, 3> delta;
	std::array<FP, 3> shift;
	for (int64_t i = 0; i < 3; i++) {
		nelem *= shape[i];
		int64_t shn1 = shape[i] - 1;
		delta[i] = shn1 == 0 ? 0.0f : size[i] / shn1;
		shift[i] = shn1 / FP(2);
	}
	Vec3<FP>* data = new Vec3<FP>[nelem];


	/* Fill arrays */
	Vec3<FP>* dptr = data;
	for (int64_t x = 0; x < shape[0]; x++) {
		for (int64_t y = 0; y < shape[1]; y++) {
			for (int64_t z = 0; z < shape[2]; z++) {
				Vec3<FP> val{ (x - shift[0]) * delta[0], (y - shift[1]) * delta[1], (z - shift[2]) * delta[2] };
				memcpy((void*)dptr, (void*)&val, sizeof(FP) * 3);
				dptr++;
			}
		}
	}
	return data;
}
template Vec3<float>* grid3d_struct_cpy(std::array<std::int64_t, 3> shape, std::array<float, 3> size);

template <typename FP>
FP* grid3d_ins(std::array<std::int64_t, 3> shape, std::array<FP, 3> size) {
	/* Half perf of og. */
	using namespace std;
	/* Initialization. */
	int64_t nelem = 1;
	std::array<FP, 3> delta;
	std::array<FP, 3> shift;
	for (int64_t i = 0; i < 3; i++) {
		nelem *= shape[i];
		int64_t shn1 = shape[i] - 1;
		delta[i] = shn1 == 0 ? 0.0f : size[i] / shn1;
		shift[i] = shn1 / FP(2);
	}
	FP* data = new FP[nelem * 3];


	/* Fill arrays */
	int64_t cpx = shape[2] * shape[1] * 3;
	for (int64_t x = 0; x < shape[0]; x++) {
		int64_t index = cpx * x;
		for (int64_t y = 0; y < shape[1]; y++) {
			for (int64_t z = 0; z < shape[2]; z++) {
				data[index++] = (x - shift[0]) * delta[0];
				data[index++] = (y - shift[1]) * delta[1];
				data[index++] = (z - shift[2]) * delta[2];
			}
		}
	}
	return data;
}
template float* grid3d_ins(std::array<std::int64_t, 3> shape, std::array<float, 3> size);


template <typename FP>
FP* grid3d_r(std::array<std::int64_t, 3> shape, std::array<FP, 3> size) {
	/* IndexRange(), no perf. difference. */
	using namespace std;
	/* Initialization. */
	int64_t nelem = 1;
	std::array<FP, 3> delta;
	std::array<FP, 3> shift;
	for (int64_t i = 0; i < 3; i++) {
		nelem *= shape[i];
		int64_t shn1 = shape[i] - 1;
		delta[i] = shn1 == 0 ? FP(0) : size[i] / shn1;
		shift[i] = shn1 / FP(2);
	}
	FP* data = new FP[nelem * 3];


	/* Fill arrays */
	FP* dptr = data;
	for (const int64_t x : IndexRange(shape[0])) {
		for (const int64_t y : IndexRange(shape[1])) {
			for (const int64_t z : IndexRange(shape[2])) {
				*dptr++ = (x - shift[0]) * delta[0];
				*dptr++ = (y - shift[1]) * delta[1];
				*dptr++ = (z - shift[2]) * delta[2];
			}
		}
	}
	return data;
}
template float* grid3d_r(std::array<std::int64_t, 3> shape, std::array<float, 3> size);


template <typename FP>
FP* grid3d_b(std::array<std::int64_t, 3> shape, std::array<FP, 3> size) {
	/* Single loop, slow */
	using namespace std;
	/* Initialization. */
	int64_t nelem = 1;
	std::array<FP, 3> delta;
	std::array<FP, 3> shift;
	for (int64_t i = 0; i < 3; i++) {
		nelem *= shape[i];
		int64_t shn1 = shape[i] - 1;
		delta[i] = shn1 == 0 ? 0.0f : size[i] / shn1;
		shift[i] = shn1 / FP(2);
	}
	FP* data = new FP[nelem * 3];

	std::array<int64_t, 3> sh_prod{ shape[0] * shape[1], shape[1], 1 };

	/* Loop over elements */
	FP* dptr = data;
	for (int64_t i = 0; i < nelem; i++) {
		std::lldiv_t dv = std::div(i, sh_prod[0]);
		*dptr++ = (dv.quot - shift[0]) * delta[0];
		dv = std::div(dv.rem, sh_prod[1]);
		*dptr++ = (dv.quot - shift[1]) * delta[1];
		*dptr++ = (dv.rem - shift[2]) * delta[2];
	}
	return data;
}
template float* grid3d_b(std::array<std::int64_t, 3> shape, std::array<float, 3> size);