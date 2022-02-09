/*
 * reformDecompose.c
 *
 *  Created on: Feb 7, 2022
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "samplingDRO.h"
configType config;

int solveReformDecompose(stocType *stoc, probType **prob, cellType *cell) {
	int candidCut;
	clock_t tic;

	cell->incumbEst = INFINITY;

	/* Main loop of the algorithm */
	while ( cell->k < config.MAX_ITER ) {
		tic = clock();
		cell->k++;

		/******* 1. Solve the subproblem with candidate solution, form and update the candidate cut *******/
		if ( (candidCut = formReformCuts(prob[1], cell, cell->candidX)) < 0 ) {
			errMsg("algorithm", "solveCell", "failed to add candidate cut", 0);
			return 1;
		}

		/******* 2. Check for optimality of the current solution *******/
		if (optimal(prob[0], cell)) {
			cell->incumbEst = cell->candidEst;
			break;
		}

		/******* 3. Solve the master problem to obtain the new candidate solution *******/
		if ( solveMaster(prob[0]->num, prob[0]->dBar, cell, prob[0]->lb) ) {
			errMsg("algorithm", "solveCell", "failed to solve master problem", 0);
			return 1;
		}

		/* Update the solution times */
		cell->time->masterAccumTime += cell->time->masterIter; cell->time->subprobAccumTime += cell->time->subprobIter;
		cell->time->argmaxAccumTime += cell->time->argmaxIter; cell->time->optTestAccumTime += cell->time->optTestIter;
		cell->time->masterIter = cell->time->subprobIter = cell->time->optTestIter = 0.0;
		cell->time->argmaxIter = cell->time->optTestIter = 0.0;
		cell->time->iterTime = ((double) clock() - tic)/CLOCKS_PER_SEC; cell->time->iterAccumTime += cell->time->iterTime;
	}//END main loop

	return 0;
}//END solveReformDecompose()

int formReformCuts(probType *prob, cellType *cell, dVector Xvect) {
	dVector piS, beta, piCBar;
	iVector indices;
	sparseMatrix COmega;
	sparseVector bOmega;
	char cutName[NAMESIZE];
	double mubBar, alpha;
	static int cummCutNum = 0;

	/* Allocate memory and initializations */
	piS = (dVector) arr_alloc(prob->num->rows+1, double);
	indices = (iVector) arr_alloc(prob->num->prevCols+2, int);

	bOmega.cnt = prob->num->rvbOmCnt; bOmega.col = prob->coord->rvbOmRows;
	COmega.cnt = prob->num->rvCOmCnt; COmega.col = prob->coord->rvCOmCols; COmega.row = prob->coord->rvCOmRows;

	for ( int i = 1; i <= prob->num->prevCols; i++ )
		indices[i] = i-1;
	indices[prob->num->prevCols+1] = prob->num->prevCols;

	/* All the subproblems are solved in any iteration. */
	/* Loops on \xi */
	for ( int m = 0; m < cell->omega->numObs; m++ ) {
		/* 1a. Construct the subproblem with a given observation and master solution, solve the subproblem, and obtain dual information. */
		if ( solveSubprob(prob, cell->subprob, Xvect, cell->omega->vals[m], &cell->spFeasFlag, &cell->time->subprobIter, piS, &mubBar) ) {
			errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
			goto TERMINATE;;
		}
		cell->LPcnt++;
		if ( ! cell->spFeasFlag ) {
			printf("Subproblem is infeasible, adding feasibility cut to the master.\n");
			break;
		}

		/* 1b. Compute the cut coefficients */
		bOmega.val = cell->omega->vals[m] + prob->coord->rvOffset[0];
		COmega.val = cell->omega->vals[m] + prob->coord->rvOffset[1];

		alpha = vXvSparse(piS, prob->bBar) + mubBar + vXvSparse(piS, &bOmega);

		beta = vxMSparse(piS, prob->Cbar, prob->num->prevCols+1);
		piCBar = vxMSparse(piS, &COmega, prob->num->prevCols);
		for (int c = 1; c <= prob->num->rvCOmCnt; c++)
			beta[prob->coord->rvCOmCols[c]] += piCBar[c];
		mem_free(piCBar);

#if defined(CUT_CHECK)
		printf("Objective function value vs. estimate = (%lf v. %lf)\n", getObjective(cell->subprob->lp, PROB_LP),
				alpha - vXv(beta, Xvect, NULL, prob->num->prevCols));
#endif

		beta[0] = 1.0;			/* coefficient of \eta */
		for ( int n = 0; n < cell->omega->numObs; n++ ) {
			/* Coefficient of \lambda */
			indices[0] = prob->num->prevCols + 1 + n;
			beta[prob->num->prevCols+1] = pNorm(cell->omega->vals[m], cell->omega->vals[n], prob->num->numRV, config.DRO_PARAM_1);

			/* Add cut to the master problem  */
			sprintf(cutName, "cut_%04d", cummCutNum++);
			/* Add the row in the solver */
			if ( addRow(cell->master->lp, prob->num->prevCols + 2, alpha, GE, 0, indices, beta, cutName) ) {
				errMsg("solver", "addcut2Master", "failed to add new row to problem in solver", 0);
				return -1;
			}
		}
		mem_free(beta);
	}//END of \xi loop

	if ( piS ) mem_free(piS);
	if ( indices ) mem_free(indices);
	return 0;

	TERMINATE:
	if ( piS ) mem_free(piS);
	if ( indices ) mem_free(indices);
	return 1;

}//END formReformCuts()
