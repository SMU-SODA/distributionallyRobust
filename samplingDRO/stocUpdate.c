/*
 * stocUpdate.c
 *
 *  Created on: Jun 24, 2018
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "samplingDRO.h"

extern configType config;

/* This function updates all the structures necessary for forming a stochastic cut. The latest observation of omega and the latest dual solution
 * to the subproblem are added to their appropriate structures. Then Pi x b and Pi x C are computed and for the latest omega and dual dVector,
 * and are added to the appropriate structures.
 * Note that the new column of delta is computed before a new row in lambda is calculated and before the new row in delta is completed,
 * so that the intersection of the new row and new column in delta is only computed once (they overlap at the bottom, right-hand corner). */
int stochasticUpdates(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, lambdaType *lambda, sigmaType *sigma,
		deltaType *delta, omegaType *omega, int iter, dVector pi, double mubBar) {
	int 	sigmaIdx, lambdaIdx;
	bool 	newLambdaFlag= false;

	/* extract the dual solutions corresponding to rows with random elements in them */
	lambdaIdx = calcLambda(num, coord, pi, lambda, &newLambdaFlag);

	/* compute Pi x bBar and Pi x Cbar */
	sigmaIdx = calcSigma(num, coord, bBar, Cbar, pi, mubBar, lambdaIdx, newLambdaFlag, iter, sigma);

	/* Only need to calculate row if a distinct lambda was found. We could use Pi, instead of lambda(Pi), for this calculation, */
	/* and save the time for expanding/reducing dVector even though the lambda is the same, the current Pi might be a
     distinct one due to the variations in sigma*/
	if (newLambdaFlag)
		calcDeltaRow(num, coord, omega, lambda, lambdaIdx, delta);

	return sigmaIdx;
}//END stochasticUpdates()

/* This function stores a new lambda_pi dVector in the lambda structure.  Each lambda_pi represents only those dual variables whose rows in the
 * constraint matrix have random elements.  Thus  the (full) dual dVector, Pi,  passed to the function is converted into the sparse dVector lambda_pi.
 * This dVector is then compared with all previous lambda_pi dVectors, searching for a duplication. If a duplicate is found, the dVector is not added
 * to the structure, and the function returns the index of the duplicate dVector. Otherwise, it adds the dVector to the end of the structure,
 *and returns an index to the last element in lambda. */
int calcLambda(numType *num, coordType *coord, dVector Pi, lambdaType *lambda, bool *newLambdaFlag) {
	int 	pi_idx;
	dVector	lambda_pi;

	/* Pull out only those elements in dual dVector which have rv's */
	lambda_pi = reduceVector(Pi, coord->rvRows, num->rvRowCnt);

	/* Compare resulting lambda_pi with all previous dVectors */
	for (pi_idx = 0; pi_idx < lambda->cnt; pi_idx++)
		if (equalVector(lambda_pi, lambda->vals[pi_idx], num->rvRowCnt, config.TOLERANCE)) {
			mem_free(lambda_pi);
			*newLambdaFlag = false;
			return pi_idx;
		}

	/* Add the dVector to lambda structure */
	lambda->vals[lambda->cnt] = lambda_pi;
	*newLambdaFlag = true;

	return lambda->cnt++;
}//END calcLambda

int calcSigma(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *CBar, dVector pi, double mubBar,
		int idxLambda, bool newLambdaFlag, int currentIter, sigmaType *sigma) {
	dVector	piCBar, temp;
	double 	pibBar;
	int 	cnt;

	/* sigma = \pi_t^\top \bar{b}_t - \bar{C}_t^\top \pi_t */
	pibBar = vXvSparse(pi, bBar) + mubBar;

	temp = vxMSparse(pi, CBar, num->prevCols);
	piCBar = reduceVector(temp, coord->CCols, num->cntCcols);
	mem_free(temp);

	if (!newLambdaFlag){
		for (cnt = 0; cnt < sigma->cnt; cnt++) {
			if (DBL_ABS(pibBar - sigma->vals[cnt].pib) <= config.TOLERANCE) {
				if (equalVector(piCBar, sigma->vals[cnt].piC, num->cntCcols, config.TOLERANCE))
					if(sigma->lambdaIdx[cnt]== idxLambda){
						mem_free(piCBar);
						return cnt;
					}
			}
		}
	}

	sigma->vals[sigma->cnt].pib  = pibBar;
	sigma->vals[sigma->cnt].piC  = piCBar;
	sigma->lambdaIdx[sigma->cnt] = idxLambda;
	sigma->ck[sigma->cnt] = currentIter;

	return sigma->cnt++;
}//END calcSigma()

/* This function calculates a new row in the delta structure, based on a new dual dVector, lambda_pi, by calculating lambda_pi X b and
 * lambda_pi X C for all previous realizations of b(omega) and C(omega).  It is assumed that the lambda dVector is distinct from all previous ones
 * and thus a new row is warranted. */
int calcDeltaRow(numType *num, coordType *coord, omegaType *omega, lambdaType *lambda, int lambdaIdx, deltaType *delta) {
	sparseVector bomega;
	sparseMatrix Comega;
	dVector 	lamb_pi, pixC;
	int		obs;

	bomega.cnt = num->rvbOmCnt;	bomega.col = coord->rvbOmRows;
	Comega.cnt = num->rvCOmCnt; Comega.col = coord->rvCOmCols; Comega.row = coord->rvCOmRows;

	if ( !(delta->vals[lambdaIdx] = (pixbCType *) arr_alloc(omega->cnt, pixbCType)))
		errMsg("allocation", "calcDeltaRow", "delta->val[cnt]", 0);

	/* expand the compressed lambda dVector */
	lamb_pi = expandVector(lambda->vals[lambdaIdx], coord->rvRows, num->rvRowCnt, num->rows);

	/* go through all the observations and compute pi x b and pi x C */
	for (obs = 0; obs < omega->cnt; obs++) {
		bomega.val= omega->vals[obs];
		Comega.val = omega->vals[obs] + num->rvbOmCnt;

		delta->vals[lambdaIdx][obs].pib = vXvSparse(lamb_pi, &bomega);
		if ( num->rvColCnt != 0 ) {
			pixC = vxMSparse(lamb_pi, &Comega, num->prevCols);
			delta->vals[lambdaIdx][obs].piC = reduceVector(pixC, coord->rvCols, num->rvColCnt);
			mem_free(pixC);
		}
		else
			delta->vals[lambdaIdx][obs].piC = NULL;
	}

	mem_free(lamb_pi);

	return 0;

}//END calcDeltaRow()

lambdaType *newLambda(int numLambda) {
	lambdaType *lambda = NULL;

	lambda = (lambdaType *) mem_malloc(sizeof(lambdaType));
	lambda->vals = (dVector *) arr_alloc(numLambda, dVector);
	lambda->cnt = 0;

	return lambda;
}//END newLambda()

sigmaType *newSigma(int numSigma) {
	sigmaType *sigma = NULL;

	sigma = (sigmaType *) mem_malloc(sizeof(sigmaType));
	sigma->lambdaIdx = (iVector) arr_alloc(numSigma, int);
	sigma->ck = (iVector) arr_alloc(numSigma, int);
	sigma->vals = (pixbCType *) arr_alloc(numSigma, pixbCType);
	sigma->cnt = 0;

	return sigma;
}//END newSigma()

deltaType *newDelta(int numDelta) {
	deltaType *delta = NULL;

	delta = (deltaType *) mem_malloc(sizeof(deltaType));
	delta->vals = (pixbCType **) arr_alloc(numDelta, pixbCType *);

	return delta;
}//END newDelta()

void freeLambdaType (lambdaType *lambda) {

	if ( lambda ) {
		if (lambda->vals ) {
			for ( int n = 0; n < lambda->cnt; n++ ) {
				if (lambda->vals[n]) mem_free(lambda->vals[n]);
			}
			mem_free(lambda->vals);
		}
		mem_free(lambda);
	}

}//END freeLambdaType()

void freeSigmaType (sigmaType *sigma) {

	if ( sigma ) {
		if ( sigma->ck ) mem_free(sigma->ck);
		if ( sigma->lambdaIdx ) mem_free(sigma->lambdaIdx);
		if ( sigma->vals ) {
			for ( int n = 0; n < sigma->cnt; n++ ) {
				if ( sigma->vals[n].piC ) mem_free(sigma->vals[n].piC);
			}
			mem_free(sigma->vals);
		}
		mem_free (sigma);
	}

}//END freeSigmaType()

void freeDeltaType (deltaType *delta, int numLambda, int numOmega) {

	if ( delta ) {
		for ( int m = 0; m < numLambda; m++ ) {
			if ( delta->vals[m] ) {
				for ( int n = 0; n < numOmega; n++ ) {
					if ( delta->vals[m][n].piC ) mem_free(delta->vals[m][n].piC);
				}
				mem_free(delta->vals[m]);
			}
		}
		mem_free(delta->vals);
		mem_free(delta);
	}

}//END freeDeltaType()

/* This function allocates memory for an omega structure.  It allocates the memory to structure elements: a dVector to hold an array of
 * observation and the probability associated with it. */
omegaType *newOmega(stocType *stoc) {
	omegaType *omega;
	int cnt, i, base, idx;

	if ( !(omega = (omegaType *) mem_malloc(sizeof(omegaType))) )
		errMsg("allocation","newOmega", "omega", 0);
	if ( !(omega->probs = (dVector) arr_alloc(config.MAX_OBS, double)) )
		errMsg("allocation", "newOmega", "omega->probs", 0);
	if ( !(omega->vals = (dVector *) arr_alloc(config.MAX_OBS, dVector)) )
		errMsg("allocation", "newOmega", "omega->vals", 0);
	omega->cnt = 0; omega->numRV = stoc->numOmega;

	if ( config.SAMPLING_TYPE == 2 ) {
		return omega;
	}
	else if ( config.SAMPLING_TYPE == 1) {
		omega->cnt = config.MAX_OBS;
		return omega;
	}

	if ( strstr(stoc->type, "BLOCKS") != NULL ) {
		if ( (omega->cnt = stoc->numVals[0]) <= config.MAX_OBS) {
			omega->vals = (dVector *) mem_realloc(omega->vals, omega->cnt*sizeof(dVector));
			omega->probs = (dVector) mem_realloc(omega->probs, omega->cnt*sizeof(double));
			for ( cnt = 0; cnt < omega->cnt; cnt++) {
				omega->probs[cnt]= stoc->probs[0][cnt];
				if ( !(omega->vals[cnt] = (dVector) arr_alloc(omega->numRV+1, double)) )
					errMsg("allocation", "updateOmega", "omega->vals[cnt]", 0);
				for (i = 0; i < omega->numRV; i++)
					omega->vals[cnt][i+1]=stoc->vals[i][cnt]-stoc->mean[i];
				omega->vals[cnt][0] = oneNorm(omega->vals[cnt]+1, omega->numRV);
			}
		}
		else {
			omega->cnt = config.MAX_OBS;
			config.SAMPLING_TYPE = 1;
		}
	}
	else if ( strstr(stoc->type, "INDEP") != NULL ) {
		omega->cnt = 1; i = 0;
		while ( i < stoc->numOmega ) {
			omega->cnt *= stoc->numVals[i];
			if (omega->cnt > config.MAX_OBS) {
				omega->cnt = config.MAX_OBS;
				config.SAMPLING_TYPE = 1;
				break;
			}
			i++;
		}

		if ( config.SAMPLING_TYPE == 0 ){
			omega->vals = (dVector *) mem_realloc(omega->vals, omega->cnt*sizeof(dVector));
			omega->probs = (dVector) mem_realloc(omega->probs, omega->cnt*sizeof(double));
			for ( cnt = 0; cnt < omega->cnt; cnt++) {
				if ( !(omega->vals[cnt] = (dVector) arr_alloc(omega->numRV+1, double)) )
					errMsg("allocation", "updateOmega", "omega->vals[cnt]", 0);
				omega->probs[cnt] = 1; base = omega->cnt;
				for ( i = 0; i < omega->numRV; i++ ) {
					base /= stoc->numVals[i];
					idx = (int)((double) cnt / (double) base) % stoc->numVals[i];
					omega->vals[cnt][i+1] = stoc->vals[i][idx]-stoc->mean[i];
					omega->probs[cnt] *= stoc->probs[i][idx];
				}
			}
		}
	}
	else {
		omega->cnt = config.MAX_OBS;
		config.SAMPLING_TYPE = 1;
	}

	return omega;
}//END newOmega()

void freeOmegaType(omegaType *omega, bool partial) {
	int n;

	if ( omega->vals ) {
		for ( n = 0; n < omega->cnt; n++ )
			if ( omega->vals[n] )
				mem_free(omega->vals[n]);
		if ( partial ) {
			omega->cnt = 0;
			return;
		}
		mem_free(omega->vals);
	}
	if ( omega->probs ) mem_free(omega->probs);
	mem_free(omega);

}//END freeOmegaType()
