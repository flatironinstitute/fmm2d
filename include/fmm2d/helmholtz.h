#pragma once

#include <complex.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

void hfmm2d_s_c_p(double eps, double _Complex zk, int64_t ns, double const *sources, double _Complex const *charge, double _Complex *pot, int64_t *ier);

#ifdef __cplusplus
}
#endif
