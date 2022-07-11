#ifndef HYPER_LOG_LOG_H_INCLUDED
#define HYPER_LOG_LOG_H_INCLUDED

#include <stdint.h>

/*
 * HyperLogLog supports sparse and dense represenations.
 * Dense representation is a classical representation:
 * there is alsways allocated (1 << precision) number of
 * counters even if some of them are not used. It should
 * be used for esimating large cardinalities.
 * Sparse representation allocates only pairs of index and its the
 * highest rank that has added. It requires less memory for
 * small cardinalities and can provide better accuracy. Sparse
 * representation swithces to dence representation if it statrs
 * to require more amount of memory that is needed for dense
 * representation.
 */
enum HLL_REPRESENTATION {
	HLL_SPARSE,
	HLL_DENSE
};

/*
 * Estimator that is used for the HyperLogLog algorithm.
 * The algorithm allows to estimate cardinality of a multiset
 * using fixed aount of memory. Memory requirements and estimation
 * accuracy are determined by the algorithm precision parameter.
 * The relative error is 1.04/sqrt(m) and the memory capasity is m*6 bits
 * where m is number of counters wich equals to 2^precision.
 */
struct hll {
	/* See the comment to HLL_REPRESENTATION enum. */
	enum HLL_REPRESENTATION representation;
	/*
	 * Interpretation of this data depends on representation.
	 * For dense representaion:
	 * Array of registers of size 2^precision.
	 * Every register stores the maximum received rank of set of
	 * hashes wich last precision bits are equal to register index.
	 * For sparse representation:
	 * Stucture of sparsely represented HyperLogLog.
	 * 
	 */
	uint8_t *data;
	/*
	 * Precision is equal to number of bits that are interpreted
	 * as register index. Available values are from HLL_MIN_PRECISION to
	 * HLL_MAX_PRECISION (defined in hll_emprirical.h.
	 * The larger value leads to less estimation error
	 * but larger memory requirement (2^precision * 6 bits).
	 */
	uint8_t precision;
	/*
	 * Cached value of the last estimation.
	 */
	double cached_estimation;
};

/*
 * Creates a HyperLogLog estimator. Precision defines
 * the estimation error and memory requirements.
 * The algorithm needs 2^precision * 6 bits memory.
 * Set precision as 14 for an estimation error of less than 1%.
 * The precision can take any value from HLL_MIN_PRECISION to
 * HLL_MAX_PRECISION (defined in hll_emprirical.h).
 */
struct hll *
hll_create(uint8_t precision);

/*
 * Add a hash of a dataset element to the hll estimator.
 */
void
hll_add(struct hll *hll, uint64_t hash);

/*
 * Estimate cardinality of the hll estimator.
 */
uint64_t
hll_estimate(struct hll *hll);

/*
 * Destroy the hll structure.
 */
void
hll_destroy(struct hll *hll);

#endif /* HYPER_LOG_LOG_H_INCLUDED */
