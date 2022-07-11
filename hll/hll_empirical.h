#ifndef HLL_EMPIRICAL_H_INCLUDED
#define HLL_EMPIRICAL_H_INCLUDED

#include <stdint.h>

/*
 * HyperLogLog needs a bias correction for small cardianlities.
 * The bias corrections was found empirically for precsion range
 * defined by following constants.
 */
enum HLL_EMPERICAL_LIMITS {
	HLL_MIN_PRECISION = 6,
	HLL_MAX_PRECISION = 18
};

/*
 * Return empirical based interpolated bias correction for the raw_estimation.
 */
double hll_empirical_bias_correction(uint8_t precision, double raw_estimation);

/*
 * Return the threshold below which linear counting algorithm has a smaller
 * error than HyperLogLog.
 */
uint64_t hll_empirical_estimation_threshold(uint8_t precision);

#endif /* HLL_EMPIRICAL_H_INCLUDED */
