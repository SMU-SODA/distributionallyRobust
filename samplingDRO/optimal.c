/*
 * optimal.c
 *
 *  Created on: Sep 15, 2020
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "samplingDRO.h"

extern configType config;

bool optimal(probType *prob, cellType *cell) {

	double est = vXvSparse(cell->candidX, prob->dBar) + maxCutHeight(cell->cuts, cell->candidX, prob->num->cols, 0);

#if defined(ALGO_CHECK)
	printf("Current estimate = %lf; Previous estimate = %lf\n", est, cell->candidEst);
#endif

	clock_t tic = clock();
	if ( config.SAMPLING_TYPE == 1 ) {
		if ( cell->k > config.MIN_ITER ) {
			if (cell->candidEst >= 0){
				cell->optFlag = (est < (1 + config.EPSILON) * cell->candidEst);
			}
			else
				cell->optFlag = (est > (1 + config.EPSILON) * cell->candidEst);
		}
	}
	else {
		fprintf(stderr, "Optimality rules not set.\n");
	}

	cell->time->optTestIter += ((double) (clock() - tic))/CLOCKS_PER_SEC;

	return cell->optFlag;
}//END optimal()

