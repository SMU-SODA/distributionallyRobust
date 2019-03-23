/*
 * algo.c
 *
 *  Created on: Mar 22, 2019
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */
#include "samplingDRO.h"

extern configType config;

int algo(oneProblem *orig, timeType *tim, stocType *stoc, cString inputDir, cString probName) {
	probType **prob = NULL;
	cellType *cell = NULL;
	dVector  meanSol;
	clock_t	 tic;
	int		 rep;

	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell, &meanSol) )
		goto TERMINATE;

	for ( rep = 0; rep < config.NUM_REPS; rep++ ) {
		fprintf(stdout, "\n====================================================================================================================================\n");
		fprintf(stdout, "Replication-%d\n", rep+1);

		/* setup the seed to be used in the current iteration */
		config.RUN_SEED[0] = config.RUN_SEED[rep+1];

		if ( rep != 0 ) {
			/* clean up the cell for the next replication */
			if ( cleanCellType(cell, prob[0], meanSol) ) {
				errMsg("algorithm", "benders", "failed to solve the cells using MASP algorithm", 0);
				goto TERMINATE;
			}
		}

		/* Update omega structure */
		if ( config.SAMPLING_TYPE == 1 ) {
			cell->omega->cnt = config.MAX_OBS;
			setupSAA(stoc, &config.RUN_SEED[0], &cell->omega->vals, &cell->omega->probs, &cell->omega->cnt, config.TOLERANCE);
			for ( int m = 0; m < cell->omega->cnt; m++ )
				for ( int n = 1; n <= stoc->numOmega; n++ )
					cell->omega->vals[m][n] -= stoc->mean[n-1];
		}

		tic = clock();
		/* Use two-stage algorithm to solve the problem */
		if ( solveFixedDROCell(stoc, prob, cell) ) {
			errMsg("algorithm", "benders", "failed to solve the cells using MASP algorithm", 0);
			goto TERMINATE;
		}
		cell->time->repTime = ((double) (clock() - tic))/CLOCKS_PER_SEC;
	}

	return 0;

	TERMINATE:
	return 1;
}//END algo()

int solveFixedDROCell(stocType *stoc, probType **prob, cellType *cell) {

	return 0;
}//END solveFixedDROCell()
