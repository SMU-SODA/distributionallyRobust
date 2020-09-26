/*
 * cuts.c
 *
 *  Created on: Mar 23, 2019
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "samplingDRO.h"

extern configType config;

int formDeterministicCut(probType *prob, cellType *cell, dVector Xvect) {
	oneCut *cut;
	dVector piS, alpha, *beta, piCBar, spObj;
	sparseMatrix COmega;
	sparseVector bOmega;
	double mubBar;
	int obs, cutID = -1;
	clock_t tic;

	/* Allocate memory and initializations */
	piS = (dVector) arr_alloc(prob->num->rows+1, double);
	spObj = (dVector) arr_alloc(cell->omega->cnt, double);
	alpha = (dVector) arr_alloc(cell->omega->cnt, double);
	beta = (dVector *) arr_alloc(cell->omega->cnt, double);

	bOmega.cnt = prob->num->rvbOmCnt; bOmega.col = prob->coord->rvbOmRows;
	COmega.cnt = prob->num->rvCOmCnt; COmega.col = prob->coord->rvCOmCols; COmega.row = prob->coord->rvCOmRows;
	cut = newCut(prob->num->prevCols, cell->k, cell->omega->cnt);

	/* All the subproblems are solved in any iteration. */
	for ( obs = 0; obs < cell->omega->cnt; obs++ ) {
		/* 1a. Construct the subproblem with a given observation and master solution, solve the subproblem, and obtain dual information. */
		if ( solveSubprob(prob, cell->subprob, Xvect, cell->omega->vals[obs], &cell->spFeasFlag, &cell->time->subprobIter, piS, &mubBar) ) {
			errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
			goto TERMINATE;;
		}
		cell->LPcnt++;

		/* Extract the objective function value of the subproblem. */
		spObj[obs] = getObjective(cell->subprob->lp, cell->subprob->type);

		if ( ! cell->spFeasFlag ) {
			printf("Subproblem is infeasible, adding feasibility cut to the master.\n");
			break;
		}

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

	tic = clock();
	/* 2. Solve the distribution separation problem if we are not solving the risk-neutral version. */
	if ( obtainProbDist(cell->sep, prob->mean, cell->omega, spObj, false, 0) ) {
		errMsg("algorithm", "formOptCut", "failed to solve the distribution separation problem", 0);
		return 1;
	}
	cell->time->distSepTime += ((double) (clock() - tic))/CLOCKS_PER_SEC;

	/* 3. Compute the aggregated cut coefficients */
	for ( obs = 0; obs < cell->omega->cnt; obs++ ) {
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

	TERMINATE:
	if ( piS ) mem_free(piS);
	if ( spObj ) mem_free(spObj);
	if ( alpha) mem_free(alpha);
	if ( beta ) {
		for ( obs = 0; obs < cell->omega->cnt; obs++ ) {
			if (beta[obs]) mem_free(beta[obs]);
		}
		mem_free(beta);
	}
	return 0;
}//END formOptCut()

int formStochasticCut(probType *prob, cellType *cell, dVector Xvect, int obsStar, bool newOmegaFlag) {
	oneCut 	*cut;
	dVector piS, spObj, piCbarX;
	double	mubBar;
	int    	cutIdx;
	clock_t tic;

	/* allocate memory to hold a subproblem duals, estimated objective function and the new cut */
	piS = (dVector) arr_alloc(prob->num->rows+1, double);
	spObj = (dVector) arr_alloc(cell->omega->cnt, double);
	cut = newCut(prob->num->prevCols, cell->omega->cnt, cell->k);

	/*********************************** 1. Subproblem solve and update resource function approximation **************************************/
	/* (a) Construct the subproblem with input observation and master solution, solve the subproblem, and complete stochastic updates */
	if ( solveSubprob(prob, cell->subprob, Xvect, cell->omega->vals[obsStar], &cell->spFeasFlag, &cell->time->subprobIter, piS, &mubBar) ) {
		errMsg("algorithm", "solveAgents", "failed to solve the subproblem", 0);
		goto TERMINATE;
	}

	if ( ! cell->spFeasFlag ) {
		printf("Subproblem is infeasible, adding feasibility cut to the master.\n");
		goto TERMINATE;
	}

	/* Increment the number of subproblems solved during algorithm */
	cell->LPcnt++;
	spObj[obsStar] = getObjective(cell->subprob->lp, cell->subprob->type);

	/* (b) Update the stochastic elements in the problem */
	tic = clock();
	cut->iStar[obsStar] = stochasticUpdates(prob->num, prob->coord, prob->bBar, prob->Cbar, cell->lambda, cell->sigma, cell->delta, config.MAX_ITER,
			cell->omega, newOmegaFlag, obsStar, cell->k, piS, mubBar);
	cell->time->argmaxIter += ((double) (clock() - tic))/CLOCKS_PER_SEC;

#ifdef STOCH_CHECK
	obj = sigma->vals[status].pib - vXv(sigma->vals[status].piC, Xvect, prob->coord->CCols, prob->num->cntCcols);
	obj += delta->vals[sigma->lambdaIdx[status]][omegaIdx].pib - vXv(delta->vals[sigma->lambdaIdx[status]][omegaIdx].piC,
			omega->vals[omegaIdx], prob->coord->rvCOmCols, prob->num->rvCOmCnt);
	printf("Objective function estimate    = %lf\n", obj);
#endif

	/* (c) Obtain the best dual vertex using the argmax procedure and the corresponding objective function estimate. */
	piCbarX = (dVector) arr_alloc(cell->sigma->cnt, double);
	double argmax;
	int istar;
	/* Pre-compute pi x Cbar x x as it is independent of observations */
	for (int c = 0; c < cell->sigma->cnt; c++)
		piCbarX[c] = vXv(cell->sigma->vals[c].piC, Xvect, prob->coord->CCols, prob->num->cntCcols);

	/* Loop through the observations to perform the argmax procedure. */
	for (int obs = 0; obs < cell->omega->cnt; obs++) {
		/* identify the maximal Pi for each observation */
		istar = computeIstar(prob->num, prob->coord, cell->sigma, cell->delta, piCbarX, Xvect,
				obs, &argmax);
		if (istar < 0) {
			errMsg("algorithm", "SDCut", "failed to identify maximal Pi for an observation", 0);
			goto TERMINATE;
		}
		cut->iStar[obs] = istar;
		spObj[obs] = argmax;
	}

	/*********************************** 2 Obtain the maximal probability distribution **************************************/
	/* Solve the distribution separation problem if we are not solving the risk-neutral version. */
	if ( config.DRO_TYPE == RISK_NEUTRAL && cell->k == 1)  {
		for ( int obs = 0; obs < cell->omega->cnt; obs++ ) {
			cell->omega->probs[obs] = cell->omega->weights[obs]/(double) cell->omega->cnt;
		}
	}
	else  {
		clock_t tic = clock();
		if ( obtainProbDist(cell->sep, prob->mean, cell->omega, spObj, obsStar, newOmegaFlag) ) {
			errMsg("algorithm", "formOptCut", "failed to solve the distribution separation problem", 0);
			return 1;
		}
		cell->time->distSepTime += ((double) (clock() - tic))/CLOCKS_PER_SEC;
	}

#if defined(SEP_CHECK)
	double objEst = 0.0, objEst_dro = 0.0;
	for (int obs = 0; obs < cell->omega->cnt; obs++) {
		objEst += cell->omega->weights[obs]*spObj[obs];
		objEst_dro += cell->omega->probs[obs]*spObj[obs];
	}
	objEst /= (double) cell->omega->numObs;
	printf("\tRisk neutral estimate = %lf, DRO estimate = %lf\n", objEst, objEst_dro);
#endif

	/*********************************** 3. Create the affine function using the stochastic elements **************************************/
	/* (a) Create an affine lower bound */
	cut->alpha = 0.0;
	cut->beta[0] = 1.0;
	for ( int obs = 0; obs < cell->omega->cnt; obs++ ) {
		cut->alpha += cell->sigma->vals[cut->iStar[obs]].pib * cell->omega->probs[obs];
		cut->alpha += cell->delta->vals[cell->sigma->lambdaIdx[cut->iStar[obs]]][obs].pib * cell->omega->probs[obs];

		for (int c = 1; c <= prob->num->cntCcols; c++)
			cut->beta[prob->coord->CCols[c]] += cell->sigma->vals[cut->iStar[obs]].piC[c] * cell->omega->probs[obs];
		for (int c = 1; c <= prob->num->rvCOmCnt; c++)
			cut->beta[prob->coord->rvCols[c]] += cell->delta->vals[cell->sigma->lambdaIdx[cut->iStar[obs]]][obs].piC[c] * cell->omega->probs[obs];
	}
	cut->numObs = cell->omega->numObs;

#if defined(CUT_CHECK)
	double est1 = 0.0;
	for (int obs = 0; obs < cell->omega->cnt; obs++) {
		est1 += cell->omega->probs[obs]*spObj[obs];
	}
	double est2 = cutHeight(cut, cell->k, Xvect, prob->num->prevCols, 0, false);
	printf("\tMaximal expectation = %lf, estimated cut height = %lf\n", est1, est2);
#endif

	/* (b) Add cut to the structure and master problem  */
	if ( (cutIdx = addCut2Master(cell, cell->cuts, cut, prob->num->prevCols, 0)) < 0 ) {
		errMsg("algorithm", "formSDCut", "failed to add the new cut to master problem", 0);
		goto TERMINATE;
	}

	mem_free(piS);
	mem_free(spObj);
	mem_free(piCbarX);

	return cutIdx;
	TERMINATE:
	return -1;
}//END formStochasticCut()

cutsType *newCuts(int maxCuts) {
	cutsType *cuts;

	if (maxCuts == 0)
		return NULL;

	if (!(cuts = (cutsType *) mem_malloc (sizeof(cutsType))))
		errMsg("allocation", "newCuts", "cuts",0);
	if (!(cuts->vals = (oneCut **) arr_alloc (maxCuts, oneCut)))
		errMsg("allocation", "newCuts", "oneCuts",0);
	cuts->cnt = 0;

	return cuts;
}//END newCuts

oneCut *newCut(int numX, int currentIter, int numSamples) {
	oneCut *cut;

	cut = (oneCut *) mem_malloc (sizeof(oneCut));
	cut->ck = currentIter;
	cut->alpha = 0.0;
	cut->alphaIncumb = 0.0;
	if (!(cut->beta = arr_alloc(numX + 1, double)))
		errMsg("allocation", "new_cut", "beta", 0);
	cut->name = (cString) arr_alloc(NAMESIZE, char);

	cut->iStar = (iVector) arr_alloc(numSamples, int);

	cut->rowNum = -1;
	cut->isIncumb = false; /* new cut is by default not an incumbent */

	return cut;
}//END newCut

/* This function loops through a set of cuts and find the highest cut height at the specified position x */
double maxCutHeight(cutsType *cuts, int currIter, dVector xk, int betaLen, double lb, bool scale) {
	double Sm = -INF, ht = 0.0;
	int cnt;

	for (cnt = 0; cnt < cuts->cnt; cnt++) {
		ht = cutHeight(cuts->vals[cnt], currIter, xk, betaLen, lb, scale);
		if (Sm < ht) {
			Sm = ht;
		}
	}

	return Sm;
}//END maxCutHeight

/* This function calculates and returns the height of a given cut at a given X.  It includes the k/(k-1) update, but does not include
 * the coefficients due to the cell. */
double cutHeight(oneCut *cut, int currObs, dVector xk, int betaLen, double lb, bool scale) {
	double height;
	double t_over_k = ((double) cut->numObs / (double) currObs);

	/* A cut is calculated as alpha - beta x X */
	height = cut->alpha - vXv(cut->beta, xk, NULL, betaLen);

	if ( scale ) {
		/* Weight cut based on number of observations used to form it */
		height *= t_over_k;

		/* Updated for optimality cut height*/
		height += (1 - t_over_k) * lb;
	}

	return height;
}//END cutHeight()

/*This function loops through all the dual vectors found so far and returns the index of the one which satisfies the expression:
 * 				argmax { Pi x (R - T x X) | all Pi }
 * where X, R, and T are given.  It is calculated in this form:
 * 				Pi x bBar + Pi x bomega + (Pi x Cbar) x X + (Pi x Comega) x X.
 * Since the Pi's are stored in two different structures (sigma and delta), the index to the maximizing Pi is actually a structure
 * containing two indices.  (While both indices point to pieces of the dual vectors, sigma and delta may not be in sync with one
 * another due to elimination of non-distinct or redundant vectors. */
int computeIstar(numType *num, coordType *coord, sigmaType *sigma, deltaType *delta, dVector piCbarX, dVector Xvect,
		int obs, double *argmax) {
	double 	arg;
	int 	cnt, maxCnt, lambdaIdx;


	*argmax = -DBL_MAX; maxCnt = 0;
	/* Run through the list of basis to choose the one which provides the best lower bound */
	for ( cnt = 0; cnt < sigma->cnt; cnt++ ) {
		lambdaIdx = sigma->lambdaIdx[cnt];
		arg = 0.0;

		arg += sigma->vals[cnt].pib + delta->vals[lambdaIdx][obs].pib - piCbarX[cnt];
		arg -= vXv(delta->vals[lambdaIdx][obs].piC, Xvect, coord->rvCOmCols, num->rvCOmCnt);

		if (arg > (*argmax)) {
			*argmax = arg;
			maxCnt = cnt;
		}
	}

	if ( (*argmax == -DBL_MAX ) )
		return -1;
	else
		return maxCnt;
}//END computeIstar

int reduceCuts() {

	return 0;
}//END reduceCuts()

void freeOneCut(oneCut *cut) {

	if (cut) {
		if (cut->beta) mem_free(cut->beta);
		if (cut->name) mem_free(cut->name);
		if (cut->iStar) mem_free(cut->iStar);
		mem_free(cut);
	}
}

void freeCutsType(cutsType *cuts, bool partial) {
	int cnt;

	for (cnt = 0; cnt < cuts->cnt; cnt++)
		freeOneCut(cuts->vals[cnt]);

	if ( partial )
		cuts->cnt = 0;
	else {
		mem_free(cuts->vals);
		mem_free(cuts);
	}
}//END freeCuts
