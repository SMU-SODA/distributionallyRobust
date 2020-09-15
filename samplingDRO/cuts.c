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

int formOptCut(probType *prob, cellType *cell, dVector Xvect) {
	oneCut *cut;
	dVector piS, alpha, *beta, piCBar, spObj;
	sparseMatrix COmega;
	sparseVector bOmega;
	double mubBar;
	int obs, cutID = -1;

	/* Allocate memory and initializations */
	if ( !(piS = (dVector) arr_alloc(prob->num->rows+1, double)) )
		errMsg("allocation", "formOptCut", "piS", 0);
	if ( !(spObj = (dVector) arr_alloc(cell->omega->cnt, double)) )
		errMsg("allocation", "formOptCut", "spObj", 0);
	if ( !(alpha = (dVector) arr_alloc(cell->omega->cnt, double)) )
		errMsg("allocation", "formOptCut", "alpha", 0);
	if ( !(beta = (dVector *) arr_alloc(cell->omega->cnt, double)) )
		errMsg("allocation", "formOptCut", "beta", 0);
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

		cell->omega->probs[obs] = 0.01;
	}

	/* 2. Solve the distribution separation problem */
	if ( obtainProbDist(cell->sep, cell->omega->probs, spObj, cell->omega->cnt) ) {
		errMsg("algorithm", "formOptCut", "failed to solve the distribution separation problem", 0);
		return 1;
	}

	/* 3. Compute the aggregated cut coefficients */
	for ( obs = 0; obs < cell->omega->cnt; obs++ ) {
		cut->alpha += cell->omega->probs[obs]*alpha[obs];
		for (int c = 1; c <= prob->num->prevCols; c++)
			cut->beta[c] += cell->omega->probs[obs]*beta[obs][c];
	}

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

oneCut *newCut(int numX, int currentIter, int numObs) {
	oneCut *cut;

	cut = (oneCut *) mem_malloc (sizeof(oneCut));
	cut->isIncumb = false; 								/* new cut is by default not an incumbent */
	cut->alphaIncumb = 0.0;
	cut->rowNum = -1;
	cut->ck = currentIter;
	if (!(cut->beta = arr_alloc(numX + 1, double)))
		errMsg("allocation", "new_cut", "beta", 0);
	cut->alpha = 0.0;
	cut->name = (cString) arr_alloc(NAMESIZE, char);
	cut->omegaID = 0;
	cut->iStar = (iVector) arr_alloc(numObs, int);

	return cut;
}//END newCut

/* This function loops through a set of cuts and find the highest cut height at the specified position x */
double maxCutHeight(cutsType *cuts, dVector xk, int betaLen, int obsID) {
	double Sm = -INF, ht = 0.0;
	int cnt;

	for (cnt = 0; cnt < cuts->cnt; cnt++) {
		if ( cuts->vals[cnt]->omegaID == obsID ) {
			ht = cutHeight(cuts->vals[cnt], xk, betaLen);
			if (Sm < ht)
				Sm = ht;
		}
	}

	return Sm;
}//END maxCutHeight

/* This function calculates and returns the height of a given cut at a given X. */
double cutHeight(oneCut *cut, dVector xk, int betaLen) {
	double height;

	/* A cut is calculated as alpha - beta x X */
	height = cut->alpha - vXv(cut->beta, xk, NULL, betaLen);

	return height;
}//END cutHeight()


/* This function will remove the oldest cut whose corresponding dual variable is zero (thus, a cut which was slack in last solution). */
int reduceCut(oneProblem *master, cutsType *cuts, dVector vectX, dVector piM, int betaLen, iVector iCutIdx,
		omegaType *omega, int obsID) {
	double 	height, minHeight;
	int 	oldestCut, idx;

	oldestCut = cuts->cnt-1;

	/* identify the oldest loose cut */
	for (idx = 0; idx < cuts->cnt; idx++) {
		if ( idx == iCutIdx[obsID] || cuts->vals[idx]->rowNum < 0 || cuts->vals[idx]->omegaID != obsID )
			/* avoid dropping incumbent cut and newly added cuts */
			continue;

		if ( cuts->vals[idx]->ck < cuts->vals[oldestCut]->ck &&
				DBL_ABS(piM[cuts->vals[idx]->rowNum + 1]) <= config.TOLERANCE ) {
			oldestCut = idx;
		}
	}

	/* if the oldest loose cut is the most recently added cut, then the cut with minimum cut height will be dropped */
	if ( oldestCut == cuts->cnt-1 ) {
		minHeight = -INF; oldestCut = 0;

		for (idx = 0; idx < cuts->cnt; idx++) {
			if (idx == iCutIdx[obsID] || cuts->vals[idx]->omegaID != obsID )
				continue;

			height = cutHeight(cuts->vals[idx], vectX, betaLen);
			if (height < minHeight) {
				minHeight = height;
				oldestCut = idx;
			}
		}
	}

	/* drop the selected cut and swap the last cut into its place */
	if ( dropCut(master, cuts, oldestCut, iCutIdx, obsID) ){
		errMsg("algorithm", "reduceCuts", "failed to drop a cut", 0);
		return -1;
	}

	return oldestCut;
}//END reduceCuts()

/* This function removes a cut from both the cutType structure and the master problem constraint matrix.  In the cuts->vals array, the last
 * cut is swapped into the place of the exiting cut.  In the constraint matrix, the row is deleted, and the row numbers of all constraints
 * below it are decremented. */
int dropCut(oneProblem *master, cutsType *cuts, int cutIdx, iVector iCutIdx, int obsID) {
	int idx, deletedRow;

	deletedRow = cuts->vals[cutIdx]->rowNum;
	/* Get rid of the indexed cut on the solver */
	if (  removeRow(master->lp, deletedRow, deletedRow) ) {
		errMsg("solver", "dropCut", "failed to remove a row from master problem", 0);
		return 1;
	}
	freeOneCut(cuts->vals[cutIdx]);

	/* move the last cut to the deleted cut's position (structure) */
	cuts->vals[cutIdx] = cuts->vals[--cuts->cnt];

	/* if the swapped cut happens to be the incumbent cut, then update its index */
	if ( cutIdx >= cuts->cnt)
		printf("What is going on?\n");
	if ( cuts->vals[cutIdx]->isIncumb ) {
		iCutIdx[cuts->vals[cutIdx]->omegaID] = cutIdx;
	}

	for (idx = 0; idx < cuts->cnt; idx++) {
		if (cuts->vals[idx]->rowNum > deletedRow)
			--cuts->vals[idx]->rowNum;
	}

	/* decrease the number of rows on solver */
	master->mar--;

	return 0;
}//END dropCut()


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
