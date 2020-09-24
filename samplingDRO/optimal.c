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
	double est;
	clock_t tic = clock();

	if ( cell->k > config.MIN_ITER ) {
		if ( config.SAMPLING_TYPE == 1 ) {
			est = vXvSparse(cell->candidX, prob->dBar) + maxCutHeight(cell->cuts, cell->k, cell->candidX, prob->num->cols, prob->lb, false);
			if (cell->candidEst >= 0){
				cell->optFlag = (est <= (1 + config.EPSILON) * cell->candidEst);
			}
			else
				cell->optFlag = (est > (1 + config.EPSILON) * cell->candidEst);
		}
		else if ( config.SAMPLING_TYPE == 2 ) {
			est = cell->candidEst;
			if (cell->candidEst >= 0){
				cell->optFlag = (est >= (1 - config.EPSILON) * cell->incumbEst);
			}
			else
				cell->optFlag = (est > (1 + config.EPSILON) * cell->incumbEst);
		}
		else {
			fprintf(stderr, "Optimality rules not set.\n");
		}
	}
	cell->time->optTestIter += ((double) (clock() - tic))/CLOCKS_PER_SEC;

#if defined(ALGO_CHECK)
	printf("Current estimate = %lf; Previous estimate = %lf\n", est, cell->candidEst);
#endif

	return cell->optFlag;
}//END optimal()

/***********************************************************************\
 ** This function determines whether the "stagewise descent property" is satiisfied.  If the current approximation of f_k
 ** gives a lower difference between the candidate and incumbent x than the previous approximation gave, then the incumbent
 ** x is updated to the candidate x, and the reference to the incumbent cut is updated as well.  The function updates
 ** _incumbChg_ to TRUE if the incumbent was updated, and FALSE otherwise.
 \***********************************************************************/
int checkImprovement(probType *prob, cellType *cell, int candidCut) {
	double  candidEst;

	/* Calculate height at new candidate x with newest cut included */
	candidEst = vXvSparse(cell->candidX, prob->dBar) + maxCutHeight(cell->cuts, cell->k, cell->candidX, prob->num->cols, prob->lb, true);
	cell->incumbEst = vXvSparse(cell->incumbX, prob->dBar) + maxCutHeight(cell->cuts, cell->k, cell->incumbX, prob->num->cols, prob->lb, true);

#ifdef ALGO_CHECK
	printf("\nCandidate estimate = %lf, Incumbent estimate = %lf",candidEst, cell->incumbEst);
#endif

	/* If we see considerable improvement, then change the incumbent */
	if ((candidEst - cell->incumbEst) < (config.R1 * cell->gamma)) {
		/* when we find an improvement, then we need to replace the incumbent x with candidate x */
		if ( replaceIncumbent(cell, candidEst, prob->num->cols) ) {
			errMsg("algorithm", "checkImprovement", "failed to replace incumbent solution with candidate", 0);
			return 1;
		}
		cell->iCutIdx = candidCut;
		cell->incumbChg = true;
		printf("+"); fflush(stdout);
	}

	return 0;
}//END checkImprovement()

int replaceIncumbent(cellType *cell, double candidEst, int numCols) {

	/* replace the incumbent solution with the candidate solution */
	copyVector(cell->candidX, cell->incumbX, numCols, 1);
	cell->incumbEst = candidEst;

	/* update the candidate cut as the new incumbent cut */
	cell->iCutUpdt = cell->k;

	/* Since incumbent solution is now replaced by a candidate, we assume it is feasible now */
	cell->infeasIncumb = false;

	/* gamma needs to be reset to 0 since there's no difference between candidate and incumbent*/
	cell->gamma = 0.0;

	return 0;
}//END replaceIncumbent()
