#include "hll.h"
#include "hll_empirical.h"

#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <limits.h>
#include <math.h>
#include <assert.h>

enum HLL_CONSTANTS {
	/*
	 * 6 bits are needed to store the number of
	 * leading zeros of 64 bit hash.
	 */
	HLL_RANK_BITS = 6,

	/* The maximum value that can be stored in HLL_RANK_BITS bits. */
	HLL_RANK_MAX = (1 << HLL_RANK_BITS) - 1,

	/* Number of bits in one bucket that stores 4 registers. */
	HLL_BUCKET_BITS = 24,
};

/*
 * Get a number whose first n bits are equal to ones.
 */
static uint64_t
hll_ones(uint8_t n)
{
	return ((UINT64_C(1) << n) - 1);
}

static int
hll_is_valid_cache(const struct hll *hll)
{
	return hll->cached_estimation >= 0;
}

static void
hll_invalidate_cache(struct hll *hll)
{
	hll->cached_estimation = -1.f;
}

/*
 * The highest precision bits of the hash are interpreted as a register index.
 */
static uint64_t
hll_hash_register_idx(uint64_t hash, uint8_t precision)
{
	assert(precision < 64);
	return hash >> (64 - precision);
}

/*
 * Return the number of leading zeros of the first
 * (64 - precision) hash bits plus one.
 */
static uint8_t
hll_hash_rank(uint64_t hash, uint8_t precision)
{
	hash |= hll_ones(precision) << (64 - precision);
	uint8_t zero_count = 0;
	uint64_t bit = 0x1;
	while ((hash & bit) == 0) {
		++zero_count;
		bit <<= 1;
	}
	uint8_t rank = zero_count + 1;
	assert(rank <= HLL_RANK_MAX);
	return rank;
}

/* Calculate the number of registers for this presision. */
static uint64_t
hll_n_registers(uint8_t precision)
{
	assert(precision <= HLL_MAX_PRECISION);
	return UINT64_C(1) << precision;
}

/* Alpha constant that HyperLogLog uses in the estimation formula. */
static double
hll_alpha(uint8_t precision)
{
	return 0.7213 / (1.0 + 1.079 / hll_n_registers(precision));
}

/* Estimate the cardinality using the LinearCounting algorithm. */
static double
linear_counting(size_t counters, size_t null_counters)
{
	return counters * log((double)counters / null_counters);
}

/*
 * ================================================================
 * There statrs the section with functions for dense representation.
 * Dense representation is a classical representation: there is alsways
 * allocated 2^precision number of counters even if some of them
 * are not used. It should be used for esimating large cardinalities.
 * ================================================================
 */

/*
 * Calculate the amount of memory reqired to store
 * the registers for dense representation.
 */
static size_t
hll_dense_reqired_memory(uint8_t precision)
{
	size_t n_registers = hll_n_registers(precision);
	return n_registers * HLL_RANK_BITS / CHAR_BIT;
}

/* Create a densely represented HyperLogLog estimator. */
struct hll *
hll_dense_create(uint8_t precision)
{
	assert(precision >= HLL_MIN_PRECISION);
	assert(precision <= HLL_MAX_PRECISION);
	struct hll *hll = calloc(1, sizeof(*hll));
	hll->representation = HLL_DENSE;
	hll->precision = precision;
	hll->cached_estimation = 0;
	/* For the dense representation data interpreted as registers. */
	const size_t registers_size = hll_dense_reqired_memory(precision);
	hll->data = calloc(registers_size, 1);
	return hll;
}

/* Destroy densely represented HyperLogLog. */
void
hll_dense_destroy(struct hll *hll)
{
	free(hll->data);
	free(hll);
}

/*
 * Dense register is represented by 6 bits if dense representation is used so
 * it can go out the range of one byte but 4 registers occupy exactly 3 bytes
 * so I called it a bucket. The registers array can always be separated to
 * this buckets because its size in bits is devided by 24 if precision more
 * than 2 (2^2 * 6 bits = 24, other sizes differ by power of 2 times).
 *
 * The structure of a bucket:
 * +----------+----------+----------+----------+
 * |0 regsiter|1 register|2 register|3 register|
 * +----------+----------+----------+----------+
 * |----------6 bits * 4 = 24 bits-------------|
 */
struct reg_bucket {
	/* Pointer to the 3 byte bucket where the register is stored. */
	uint8_t *addr;
	/* Offset of the register in the bucket. */
	size_t offset;
};

/*
 * Init a register bucket structure that is used for
 * convenient work with registers.
 */
static void
reg_bucket_init(struct reg_bucket *bucket, uint8_t *regs, size_t reg_idx)
{
	/*
	 * ASCII-visualization of the logic of the function:
	 *
	 * regs		  1 byte	 2 byte	       3 byte	      4 byte
	 * |		  |		 |	       |	      |
	 * +----------+----------+----------+----------+----------+----------+--
	 * |0 regsiter|1 register|2 register|3 register|4 register|5 register|..
	 * +----------+----------+----------+----------+----------+----------+--
	 * |	      6		 12	    18	       |          30	     32
	 * 0 bucket				       1 bucket
	 * For the 5th register:
	 * bucket_idx = 5*6 / 24 = 1,
	 * register_offset = 5*6 % 24 = 6.
	 */
	size_t bucket_bytes = HLL_BUCKET_BITS / CHAR_BIT;
	size_t bucket_idx = reg_idx * HLL_RANK_BITS / HLL_BUCKET_BITS;
	bucket->addr = (uint8_t *)(regs + bucket_idx * bucket_bytes);
	bucket->offset = reg_idx * HLL_RANK_BITS % HLL_BUCKET_BITS;
	assert(bucket->offset <= HLL_BUCKET_BITS - HLL_RANK_BITS);
}

/* Get an integer value of 3 bytes stored in the bucket. */
static uint32_t
reg_bucket_value(const struct reg_bucket *bucket)
{
	uint32_t *value_addr = (uint32_t *)bucket->addr;
	uint32_t value_mask = hll_ones(HLL_BUCKET_BITS);
	return *value_addr & value_mask;
}

/*
 * Get a mask that clears the register stored in the bucket and
 * saves the other boundary registers in the bucket.
 */
static uint32_t
reg_bucket_boundary_mask(const struct reg_bucket *bucket)
{
	/*
	 * |000000000000000000111111|
	 * |------------regstr------|
	 */
	uint32_t ones = hll_ones(HLL_RANK_BITS);
	/*
	 * |000000000000111111000000|
	 * |------------regstr------|
	 */
	uint32_t register_mask = ones << bucket->offset;
	/*
	 * |111111111111000000111111|
	 * |------------regstr------|
	 */
	uint32_t boundary_mask = ~register_mask;
	return boundary_mask;
}

/* Get the value of the register with the idx index. */
static uint8_t
hll_dense_register_rank(const struct hll *hll, size_t idx)
{
	struct reg_bucket bucket;
	reg_bucket_init(&bucket, hll->data, idx);
	uint32_t reg_mask = hll_ones(HLL_RANK_BITS);
	uint32_t bucket_value = reg_bucket_value(&bucket);
	uint8_t rank = (bucket_value >> bucket.offset) & reg_mask;
	assert(rank <= HLL_RANK_MAX);
	return rank;
}

/* Set rank as the new value of the register with the idx index. */
static void
hll_dense_set_register_rank(struct hll *hll, size_t idx, uint8_t rank)
{
	struct reg_bucket bucket;
	reg_bucket_init(&bucket, hll->data, idx);
	uint32_t boundary_mask = reg_bucket_boundary_mask(&bucket);
	uint32_t bucket_value = reg_bucket_value(&bucket);
	union {
		uint32_t value;
		uint8_t bytes[3];
	} new_bucket;
	new_bucket.value = (rank << bucket.offset) |
			   (bucket_value & boundary_mask);

	bucket.addr[0] = new_bucket.bytes[0];
	bucket.addr[1] = new_bucket.bytes[1];
	bucket.addr[2] = new_bucket.bytes[2];
}

/* Add hash to dense HyperLogLog estimator. */
void
hll_dense_add(struct hll *hll, uint64_t hash)
{
	uint8_t precision = hll->precision;
	size_t idx = hll_hash_register_idx(hash, precision);
	assert(idx < hll_n_registers(precision));
	uint8_t hash_rank = hll_hash_rank(hash, precision);
	uint8_t reg_rank = hll_dense_register_rank(hll, idx);
	if (reg_rank < hash_rank) {
		hll_dense_set_register_rank(hll, idx, hash_rank);
		hll_invalidate_cache(hll);
	}
}

/*
 * Estimate the cardinality of the dense HyperLogLog using the
 * estimation formula. Raw estimation can have larger relative error
 * for small cardinalities.
 */
static double
hll_dense_raw_estimate(const struct hll *hll)
{
	double sum = 0;
	const size_t n_registers = hll_n_registers(hll->precision);
	for (size_t i = 0; i < n_registers; ++i) {
		sum += pow(2, -hll_dense_register_rank(hll, i));
	}

	const double alpha = hll_alpha(hll->precision);
	return alpha * n_registers * n_registers / sum;
}

/* Count the number of registers that are zero. */
static size_t
hll_dense_count_zero_registers(const struct hll *hll)
{
	size_t count = 0;
	const size_t n_registers = hll_n_registers(hll->precision);
	for (size_t i = 0; i < n_registers; ++i) {
		if (hll_dense_register_rank(hll, i) == 0)
			++count;
	}
	return count;
}

/* Estimate caridnality of the dense HyperLogLog */
uint64_t
hll_dense_estimate(struct hll *hll)
{
	if (hll_is_valid_cache(hll)) {
		return hll->cached_estimation;
	}
	const uint8_t prec = hll->precision;
	const size_t n_registers = hll_n_registers(prec);
	double raw_estimation = hll_dense_raw_estimate(hll);

	double hll_estimation = raw_estimation;
	if (raw_estimation < 4.f * n_registers) {
		hll_estimation -=
			hll_empirical_bias_correction(prec, raw_estimation);
	}

	size_t zero_count = hll_dense_count_zero_registers(hll);
	double lc_estimation = zero_count != 0 ?
		linear_counting(n_registers, zero_count) :
		hll_estimation;

	uint64_t threshold = hll_empirical_estimation_threshold(prec);
	size_t estimation = lc_estimation < threshold ? lc_estimation :
							hll_estimation;
	hll->cached_estimation = estimation;
	return estimation;
}
#if 0
/*
 * ==================================================================
 * There statrs the section with functions for sparse representation.
 * Sparse representation allocates only pairs of index and its the
 * highest rank that has added. It requires less memory for
 * small cardinalities than dense representation and can provide
 * better accuracy. Sparse representation swithces to dence representation
 * if it statrs to require more amount of memory that is needed for dense
 * representation.
 * ==================================================================
 */

/*
 * Instead of registers sparse representation keeps pairs of
 * register index and its the highest rank. It helps to reduse
 * amount of memory that is needed for small cardinalities.
 * Pairs for sparse representation have the following structure:
 * First 6 bits are interpreted as its rank.
 * The last 26 bits are interpreted as its index.
 */
typedef uint32_t pair_t;

static const uint8_t HLL_SPARSE_RPECISION = 26;

static pair_t
hll_sparse_new_pair(size_t idx, uint8_t rank)
{
	pair_t pair = rank;
	pair |= idx << HLL_RANK_BITS;
	return pair;
}

static uint32_t
hll_sparse_pair_idx(pair_t pair)
{
	return pair << HLL_RANK_BITS;
}

static uint8_t
hll_sparse_pair_rank(pair_t pair)
{
	return pair & hll_ones(HLL_RANK_BITS);
}

struct pair_list_header {
	/* Number of pairs stored in the list. */
	uint32_t size;
	/* Amount of memory that is used to store the list. */
	uint32_t capacity;
};

static uint32_t
hll_sparse_pair_list_size(const uint8_t *data)
{
	struct pair_list_header *header = (struct pair_list_header *)data;
	return header->size;
}

static uint32_t
hll_sparse_pair_list_capacity(const uint8_t *data)
{
	struct pair_list_header *header = (struct pair_list_header *)data;
	return header->capacity;
}

static pair_t *
hll_sparse_pair_list(uint8_t *data)
{
	pair_t *pairs = data + sizeof(struct pair_list_header);
	return pairs;
}

static size_t
hll_sparse_memory_usage(const struct hll *hll)
{
	
}

#endif

void
hll_add(struct hll *hll, uint64_t hash)
{
	hll_dense_add(hll, hash);
#if 0
	switch (hll->representation) {
		case HLL_SPARSE:
			hll_sparse_add(hll, hash);
			break;
		case HLL_DENSE:
			hll_dense_add(hll, hash);
			break;
		default:
			unreachable();
	}
#endif
}

uint64_t
hll_estimate(struct hll *hll)
{
	return hll_dense_estimate(hll);
#if 0
	switch (hll->representation) {
		case HLL_SPARSE:
			return hll_sparse_estimate(hll);
		case HLL_DENSE:
			return hll_dense_estimate(hll);
		default:
			unreachable();
	}
#endif
}

void
hll_destroy(struct hll *hll)
{
	hll_dense_destroy(hll);
#if 0
	switch (hll->representation) {
		case HLL_SPARSE:
			hll_sparse_destroy(hll);
			break;
		case HLL_DENSE:
			hll_dense_destroy(hll);
			break;
		default:
			unreachable();
	}
#endif
}

struct hll *
hll_create(uint8_t precision)
{
	return hll_dense_create(precision);
}
