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

int solveReformDecompose(stocType *stoc, probType **prob, cellType *cell);

int solveReformDecompose(stocType *stoc, probType **prob, cellType *cell) {
	int candidCut;
	clock_t tic;

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
	oneCut *cut;
	dVector piS, alpha, *beta, piCBar, spObj;
	sparseMatrix COmega;
	sparseVector bOmega;
	double mubBar;
	int obs, cutID = -1;
	clock_t tic;

	/* Allocate memory and initializations */
	piS = (dVector) arr_alloc(prob->num->rows+1, double);
	spObj = (dVector) arr_alloc(cell->omega->numObs, double);
	alpha = (dVector) arr_alloc(cell->omega->numObs, double);
	beta = (dVector *) arr_alloc(cell->omega->numObs, double);

	bOmega.cnt = prob->num->rvbOmCnt; bOmega.col = prob->coord->rvbOmRows;
	COmega.cnt = prob->num->rvCOmCnt; COmega.col = prob->coord->rvCOmCols; COmega.row = prob->coord->rvCOmRows;
	cut = newCut(prob->num->prevCols, cell->k, cell->omega->numObs);

	/* All the subproblems are solved in any iteration. */
	for ( obs = 0; obs < cell->omega->numObs; obs++ ) {
		/* 1a. Construct the subproblem with a given observation and master solution, solve the subproblem, and obtain dual information. */
		if ( solveSubprob(prob, cell->subprob, Xvect, cell->omega->vals[obs], &cell->spFeasFlag, &cell->time->subprobIter, piS, &mubBar) ) {
			errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
			goto TERMINATE;;
		}
		cell->LPcnt++;
		if ( ! cell->spFeasFlag ) {
			printf("Subproblem is infeasible, adding feasibility cut to the master.\n");
			break;
		}

		/* Extract the objective function value of the subproblem. */
		spObj[obs] = getObjective(cell->subprob->lp, cell->subprob->type);

		/* 1b. Compute the cut coefficients */
		bOmega.val = cell->omega->vals[obs] + prob->coord->rvOffset[0];
		COmega.val = cell->omega->vals[obs] + prob->coord->rvOffset[1];

		alpha[obs] = vXvSparse(piS, prob->bBar) + mubBar + vXvSparse(piS, &bOmega);

		beta[obs] = vxMSparse(piS, prob->Cbar, prob->num->prevCols);
		piCBar = vxMSparse(piS, &COmega, prob->num->prevCols);
		for (int c = 1; c <= prob->num->rvCOmCnt; c++)
			beta[obs][prob->coord->rvCOmCols[c]] += piCBar[c];
		mem_free(piCBar);

#if defined(CUT_CHECK)
		printf("Objective function value vs. estimate = (%lf v. %lf)\n", spObj[obs], alpha[obs] - vXv(beta[obs], Xvect, NULL, prob->num->prevCols));
#endif
	}

	/* 2. Compute the aggregated cut coefficients */
	for ( obs = 0; obs < cell->omega->numObs; obs++ ) {
		cut->alpha += cell->omega->probs[obs]*alpha[obs];
		for (int c = 1; c <= prob->num->prevCols; c++)
			cut->beta[c] += cell->omega->probs[obs]*beta[obs][c];
	}
	cut->numObs = cell->omega->numObs;

	/* 4. Add cut to master. */
	cut->beta[0] = 1.0;
	if ( config.MASTER_TYPE == PROB_QP )
		cut->alphaIncumb = cut->alpha - vXv(cut->beta, cell->incumbX, NULL, prob->num->prevCols);
	else
		cut->alphaIncumb = cut->alpha;

	/* Add cut to the master problem  */
	if ( (cutID = addCut2Master(cell, cell->cuts, cut, prob->num->prevCols, 0)) < 0 ) {
		errMsg("algorithm", "formSDCut", "failed to add the new cut to master problem", 0);
		goto TERMINATE;
	}

	if ( piS ) mem_free(piS);
	if ( spObj ) mem_free(spObj);
	if ( alpha) mem_free(alpha);
	if ( beta ) {
		for ( obs = 0; obs < cell->omega->numObs; obs++ ) {
			if (beta[obs]) mem_free(beta[obs]);
		}
		mem_free(beta);
	}
	return 0;
	TERMINATE:
	if ( piS ) mem_free(piS);
	if ( spObj ) mem_free(spObj);
	if ( alpha) mem_free(alpha);
	if ( beta ) {
		for ( obs = 0; obs < cell->omega->numObs; obs++ ) {
			if (beta[obs]) mem_free(beta[obs]);
		}
		mem_free(beta);
	}
	return 1;
}//END formReformCuts()
