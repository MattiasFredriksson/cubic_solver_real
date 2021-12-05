#pragma once

#include<array>
#include<cstdint>
#include<assert.h>

template <typename FP = double, std::int64_t ndim>
FP* ndgrid(std::array<std::int64_t, ndim> shape, std::array<FP, ndim> size) {
	assert(ndim > 1);
	using namespace std;
	/* Initialization. */
	int64_t nelem = 1;
	std::array<FP, ndim> delta;
	std::array<FP, ndim> shift;
	for (int64_t i = 0; i < ndim; i++) {
		nelem *= shape[i];
		int64_t shn1 = shape[i] - 1;
		delta[i] = shn1 == 0 ? 0.0f : size[i] / shn1;
		shift[i] = shn1 / FP(2);
	}
	FP* data = new FP[nelem * ndim];

	std::array<int64_t, ndim> sh_prod;
	sh_prod[ndim - 1] = 1;
	for (int64_t i = ndim - 1; i >= 1; i--) {
		sh_prod[i - 1] = sh_prod[i] * shape[i];
	}

	/* Loop over each non-last dim in one loop. */
	int64_t nloop = (nelem + sh_prod[ndim - 2] - 1) / sh_prod[ndim - 2]; /* Integer arithmetic for: ceil(nelem / sh_prod[ndim - 2])*/
	const size_t cpy_size = ndim * sizeof(FP);
	if (ndim > 2) {
#pragma omp parallel for schedule(static)
		/* TODO: Threads should spawn over an arbitrary number of dimension(s) not limited to N - 3 (in case workload for the last 2 are minimal)*/
		for (int64_t tindex = 0; tindex < shape[ndim - 3]; tindex++) {
			FP* dptr = data + (sh_prod[ndim - 3] * tindex);
			std::array<FP, ndim> values;
			std::lldiv_t dv = std::div(sh_prod[ndim - 3] * tindex, sh_prod[0]);
			values[0] = (dv.quot - shift[0]) * delta[0];
			for (int64_t i = 1; i < ndim - 2; i++) {
				dv = std::div(dv.rem, sh_prod[i]);
				values[i] = (dv.quot - shift[i]) * delta[i];
			}
			for (int64_t y = 0; y < shape[ndim - 2]; y++) {
				values[ndim - 2] = (y - shift[ndim - 2]) * delta[ndim - 2];
				/* Iterate over last dim and assign. */
				for (int64_t x = 0; x < shape[ndim - 1]; x++) {
					values[ndim - 1] = (x - shift[ndim - 1]) * delta[ndim - 1];
					std::memcpy(static_cast<void*>(dptr), static_cast<void*>(values.data()), cpy_size);
					dptr += ndim;
				}
			}
		}
	}
	else if (false) {
		for (int64_t lindex = 0; lindex < nelem; lindex += sh_prod[ndim - 2]) {
			FP* dptr = data + lindex;
			std::array<FP, ndim> values;
			std::lldiv_t dv = std::div(lindex, sh_prod[0]);
			values[0] = (dv.quot - shift[0]) * delta[0];
			for (int64_t i = 1; i < ndim - 1; i++) {
				dv = std::div(dv.rem, sh_prod[i]);
				values[i] = (dv.quot - shift[i]) * delta[i];
			}
			/* Iterate over last dim and assign. */
			for (int64_t i = 0; i < shape.back(); i++) {
				values.back() = (i - shift.back()) * delta.back();
				std::memcpy(static_cast<void*>(dptr), static_cast<void*>(values.data()), cpy_size);
				dptr += ndim;
			}
		}
	}
	else {
		FP* dptr = data;
		std::array<FP, ndim> values;
		for (int64_t lindex = 0; lindex < nelem; lindex += sh_prod[ndim - 2]) {
			std::lldiv_t dv = std::div(lindex, sh_prod[0]);
			values[0] = (dv.quot - shift[0]) * delta[0];
			for (int64_t i = 1; i < ndim - 1; i++) {
				dv = std::div(dv.rem, sh_prod[i]);
				values[i] = (dv.quot - shift[i]) * delta[i];
			}
			/* Iterate over last dim and assign. */
			for (int64_t i = 0; i < shape.back(); i++) {
				values.back() = (i - shift.back()) * delta.back();
				std::memcpy(reinterpret_cast<void*>(dptr), reinterpret_cast<void*>(values.data()), cpy_size);
				dptr += ndim;
			}
		}
	}
	return data;
}
// template float* ndgrid<float, 3>(std::array<std::int64_t, 3> shape, std::array<float, 3> size);
