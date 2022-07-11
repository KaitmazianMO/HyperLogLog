#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <assert.h>

#include "hll.h"
#include "hll_empirical.h"

double
average_sum_of(double *arr, size_t n)
{
	double sum = 0;
	for (size_t i = 0; i < n; ++i)
		sum += arr[i];
	return sum / n;
}

double
max_of(double *arr, size_t n)
{
	double max = arr[0];
	for (size_t i = 0; i < n; ++i)
		if (max < arr[i])
			max = arr[i];
	return max;
}

double
dispersion_of(double *arr, double val, size_t n)
{
	double sqr_sum = 0;
	for (size_t i = 0; i < n; ++i) {
		sqr_sum += (val - arr[i]) * (val - arr[i]);
	}
	return sqrt(sqr_sum / (n - 1));
}

#define MAYBE_PRINT(file, ...)			\
do {						\
	if (file) {				\
		fprintf(file, __VA_ARGS__);	\
	}					\
} while (0)

struct est_errors {
	double exp_err;
	double std_err;
	double max_err;
};

uint64_t
phash(uint64_t val)
{
	unsigned char *bytes = (void *)&val;
	uint64_t p = 29791;
	uint64_t pi = 29791;
	uint64_t hash = 0;
	for (int i = 0; i < 8; ++i, pi *= p) {
		hash += pi * bytes[i];
	}
	return hash;
}

uint64_t
hash(uint64_t val)
{
	return phash(phash(val));
}

uint64_t rand_uint64()
{
	/*
	 * C standard rand() genarates random numbers in
	 * the range form 0 to RAND_MAX. RAND_MAX is at
	 * least 32767. This function helps to avoid repeats
	 * if rand() was called more than RAND_MAX times.
	 */
#if RAND_MAX < (1ULL << 16)
	uint64_t r1 = rand();
	uint64_t r2 = rand();
	uint64_t r3 = rand();
	uint64_t r4 = rand();
	uint64_t r5 = rand();
	return r1 * r2 * r3 * r4 * r5;
#elif RAND_MAX < (1ULL << 32)
	uint64_t r1 = rand();
	uint64_t r2 = rand();
	uint64_t r3 = rand();
	return r1 * r2 * r3;
#else
	uint64_t r1 = rand();
	uint64_t r2 = rand();
	return r1 * r2;
#endif
}

void
measure_hll_estimation_error(int prec, size_t min_card, size_t max_card,
			     size_t n_points, size_t sets_per_point,
			     struct est_errors *res, FILE *output)
{
	double max_err_sum = 0;
	double std_err_sum = 0;
	const size_t card_step = (max_card - min_card) / n_points;
	for (size_t n = 0; n <= n_points; ++n) {

		size_t card = min_card + card_step * n;
		double error[sets_per_point];
		double est[sets_per_point];

		for (size_t i = 0; i < sets_per_point; ++i) {

			struct hll *hll = hll_create(prec);

			for (size_t j = 0; j < card; ++j) {
				uint64_t val = rand_uint64();
				hll_add(hll, hash(val));
			}
			double this_est = hll_estimate(hll);
			double this_err = abs(this_est - card);
			est[i] = this_est;
			error[i] = this_err;

			hll_destroy(hll);
		}

		double max_err = max_of(error, sets_per_point) / (card + 1);
		max_err_sum += max_err;

		double avg_est = average_sum_of(est, sets_per_point);
		double std_err = dispersion_of(est, card, sets_per_point) / (card + 1);
		std_err_sum += std_err;

		MAYBE_PRINT(output,
			"%2d, %12zu, %12.2f, %12lg, %12lg\n",
			prec, card, avg_est, std_err, max_err);
	}

	double avg_std_err = std_err_sum / n_points;
	double avg_max_err = max_err_sum / n_points;

	res->std_err = avg_std_err;
	res->max_err = avg_max_err;
	res->exp_err = 1.04f / sqrt(1u << prec);
}

#define MAX(l, r)  ((l) < (r) ? (r) : (l))

/*
 * To reduse the running time firstly ran the test with
 * a low number of points and in case of failure for
 * some precision rerun it with larger number of points.
 */
const size_t low_n_points = 30;
const size_t large_n_points = 150;

/* Number of randomly generated sets for every cardinality. */
const size_t sets_per_point = 30;

int
main(int argc, char *argv[])
{
	/*
	 * This test can dump the data that is used to measure
	 * the estimation error. These data can be used for further
	 * analysis and empirical based impovements of the algorithm.
	 * You can increase the number of points for more accurate
	 * analysis, but it can lead to much longer execution time.
	 */
	if (argc != 1 && argc != 3) {
		fprintf(stderr,
			"Usage: %s or %s <n_points> <output_data_file>\n",
				argv[0], argv[0]);
		exit(EXIT_FAILURE);
	}

	FILE *output = NULL;
	size_t start_n_points = low_n_points;
	if (argc == 3) {
		output = fopen(argv[2], "w");
		assert(output);
		start_n_points = strtoull(argv[1], NULL, 10);
		assert(errno == 0);
	}

	struct est_errors errors[HLL_MAX_PRECISION + 1];

	MAYBE_PRINT(output,
		"prec,       card,      avg_est,      std_err,      max_err\n");

	for (int prec = HLL_MIN_PRECISION;
	     prec <= HLL_MAX_PRECISION; ++prec) {
		size_t n_regs = 1u << prec;
		size_t min_card = 0;
		size_t max_card = 10 * n_regs;
		measure_hll_estimation_error(prec, min_card, max_card,
			start_n_points, sets_per_point, errors + prec, output);
	}

	/* If some errors are to large, run the test again with more points. */
	for (int prec = HLL_MIN_PRECISION;
	     prec <= HLL_MAX_PRECISION; ++prec) {
		if (errors[prec].exp_err < errors[prec].std_err) {
			size_t n_regs = 1u << prec;
			size_t min_card = 0;
			size_t max_card = 10 * n_regs;
			measure_hll_estimation_error(prec, min_card, max_card,
				large_n_points, sets_per_point, errors + prec,
				NULL);
		}
	}

	for (int prec = HLL_MIN_PRECISION;
	     prec <= HLL_MAX_PRECISION; ++prec) {
		MAYBE_PRINT(output, "prec:%d, std_err:%lg, max_err:%lg, exp_err: %lg\n",
			prec, errors[prec].std_err, errors[prec].max_err,
			errors[prec].exp_err);
	}

	for (int prec = HLL_MIN_PRECISION;
	     prec <= HLL_MAX_PRECISION; ++prec) {
		/*
		 * The error of HyperLogLog is close to 1/sqrt(n_counters),
		 * but for small cardinalities LinearCounting is used because
		 * it has better accuracy, so the resulting error must be smaller
		 * than the HyperLogLog theoretical error.
		 */
		assert(errors[prec].std_err < errors[prec].exp_err);
	}

        return 0;
}
