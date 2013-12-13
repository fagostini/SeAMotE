#define NUM_THREADS	2	/* This number MUST NOT change!!! */
#define MAX_ID 100
#define MAX_SEQ 15000
#define MIN_MOT 3
#define MAX_MOT 21
#define NUM_BOOT 10000	/* Gives the precision of the p-value */
#define NO_NEGA 1		/* Is the factor size (Neg/Pos) that will be use to generate the negative, if not given. (USE INTEGER ONLY) */  
#define NUM_RAND 1		/* Raise this number to generate X sets of randomized sequences (NOT IN USE AT THE MOMENT) */
#define NUM_SHUFFLE 1	/* Raise this number to generate X sets of shuffled sequences (NOT IN USE AT THE MOMENT) */
#define P_TEST 1		/* Raise this number to increase the number of shuffle and random permutations (NOT IN USE AT THE MOMENT) */
#define PERCENTILE 0.05 /* NOW INITIALIZED IN THE MAIN SCRIPT. Each loop selects the top % of the motif coverage distributions ( Positive and Negative, separately ) */
#define TH_STEP 0.1		/* If the motifs tested do not achieved the set coverage, the latter will be lowered by this factor */
#define SEED_NOTA 4
#define NOTATION 14
