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
		deltaType *delta, int deltaRowLength, omegaType *omega, bool newOmegaFlag, int omegaIdx, int iter, dVector pi, double mubBar) {
	int 	sigmaIdx, lambdaIdx;
	bool 	newLambdaFlag= false;

	if ( newOmegaFlag ) {
		calcDelta(num, coord, lambda, delta, deltaRowLength, omega, newOmegaFlag, omegaIdx);
	}

	/* extract the dual solutions corresponding to rows with random elements in them */
	lambdaIdx = calcLambda(num, coord, pi, lambda, &newLambdaFlag);

	/* compute Pi x bBar and Pi x Cbar */
	sigmaIdx = calcSigma(num, coord, bBar, Cbar, pi, mubBar, lambdaIdx, newLambdaFlag, iter, sigma);

	/* Only need to calculate row if a distinct lambda was found. We could use Pi, instead of lambda(Pi), for this calculation, */
	/* and save the time for expanding/reducing dVector even though the lambda is the same, the current Pi might be a
     distinct one due to the variations in sigma*/
	if (newLambdaFlag)
		calcDelta(num, coord, lambda, delta, deltaRowLength, omega, false, lambdaIdx);

	return sigmaIdx;
}//END stochasticUpdates()

/* This function obtains a new vector of realizations of the random variables. It compares the new vector with all previous vectors, looking for
 * a duplication.  If it finds a duplicate, it returns the index of that duplicate; otherwise, it adds the vector to the list of distinct realizations
 * and returns the index of that realization. Note that the simulated observation does not have contain one-norm, while the values stored in
 * omegaType do */
int calcOmega(dVector observ, dVector stocMean, omegaType *omega, bool *newOmegaFlag, int k, double TOLERANCE) {
	int cnt;

	/* Since the problem already has the mean values on the right-hand side, remove it from the original observation */
	observ[0] = 0.0;
	for ( int m = 1; m <= omega->numOmega; m++ ) {
		observ[m] -= stocMean[m];
		observ[0] += observ[m];
	}

	/* Compare vector with all the previous observations */
	for (cnt = 0; cnt < omega->numObs; cnt++) {
		if (equalVector(observ, omega->vals[cnt], omega->numOmega, TOLERANCE)) {
			break;
		}
	}

	if ( cnt == omega->numObs ) {
		/* New observation encountered, then add the realization vector to the list. */
		omega->vals[omega->numObs] = duplicVector(observ, omega->numOmega+1);
		omega->weights[omega->numObs] = 1;
		(*newOmegaFlag) = true;
		omega->numObs++;
	}
	else {
		omega->weights[cnt]++;
		(*newOmegaFlag) = false;
	}

	/* Update omega statistics: sample mean */
	for ( int rv = 1; rv <= omega->numOmega; rv++ ) {
		for ( int m = 0; m < omega->numStats; m++ ) {
			omega->sampleStats[m][rv] = ( (k-1)*omega->sampleStats[m][rv] +
					pow(stocMean[rv] + omega->vals[cnt][rv], m+1))/(double) k;
		}
	}

#ifdef STOCH_CHECK
	printf("Observation (%d): ", *newOmegaFlag);
	printVector(omega->vals[cnt], omega->numOmega, NULL);
#endif

	return cnt;
}//END calcOmega()

/* This subroutine is used to refine the elements of omega structure after sampling a batch of observations generated using either internal or
 * external simulator.
 */
void refineOmega(omegaType *omega, dVector stocMean) {

	/* Compute sample mean */
	computeSampleStats(omega->vals, omega->weights, omega->numOmega, omega->numObs, omega->numObs, omega->sampleStats, omega->numStats);

	/* Adjust observations to be centered around mean */
	for ( int m = 0; m < omega->numObs; m++ ) {
		for ( int n = 1; n <= omega->numOmega; n++ )
			omega->vals[m][n] -= stocMean[n];
	}

	return;
}//END refineOmega()

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

/* This function calculates a new column in the delta structure, based on a new observation of omega. Thus, lambda_pi X C and lambda_pi X b
 * are calculated for all values of lambda_pi, for the new C(omega) and b(omega).  Room in the array has already been allocated, so the function
 * only fills it, in the column specified by _obs_. It is assumed that this observation is distinct from all previous ones, and thus a new column
 * must be calculated. */
int calcDelta(numType *num, coordType *coord, lambdaType *lambda, deltaType *delta, int deltaRowLength, omegaType *omega, bool newOmegaFlag, int elemIdx) {
	sparseMatrix COmega;
	sparseVector bOmega;
	dVector		 lambdaPi, piCrossC;
	int 		 idx;

	/* extract the coordinates and number of random elements */
	bOmega.cnt = num->rvbOmCnt;	bOmega.col = coord->rvbOmRows;
	COmega.cnt = num->rvCOmCnt; COmega.col = coord->rvCOmCols; COmega.row = coord->rvCOmRows;

	if ( newOmegaFlag ) {
		/* Case I: New observation encountered. */
		bOmega.val= omega->vals[elemIdx];
		COmega.val = omega->vals[elemIdx] + num->rvbOmCnt;

		/* For all dual vectors, lambda(pi), calculate pi X bomega and pi X Comega */
		for (idx = 0; idx < lambda->cnt; idx++) {
			/* Retrieve a new (sparse) dual vector, and expand it into a full vector */
			lambdaPi = expandVector(lambda->vals[idx], coord->rvRows, num->rvRowCnt, num->rows);

			/* Multiply the dual vector by the observation of bomega and Comega */
			/* Reduce PIxb from its full vector form into a sparse vector */
			delta->vals[idx][elemIdx].pib = vXvSparse(lambdaPi, &bOmega);
			if ( num->rvCOmCnt != 0 ) {
				piCrossC = vxMSparse(lambdaPi, &COmega, num->prevCols);
				delta->vals[idx][elemIdx].piC = reduceVector(piCrossC, coord->rvCOmCols, num->rvCOmCnt);
				mem_free(piCrossC);
			}
			else
				delta->vals[idx][elemIdx].piC = NULL;

			mem_free(lambdaPi);
		}
	}
	else {
		/* Case II: New dual vector encountered. */
		if ( !(delta->vals[elemIdx] = (pixbCType *) arr_alloc(deltaRowLength, pixbCType)))
			errMsg("allocation", "calcDeltaRow", "delta->val[cnt]", 0);

		/* expand the compressed lambda vector */
		lambdaPi = expandVector(lambda->vals[elemIdx], coord->rvRows, num->rvRowCnt, num->rows);

		/* go through all the observations and compute pi x b and pi x C */
		for (idx = 0; idx < omega->numObs; idx++) {

			bOmega.val = omega->vals[idx];
			COmega.val = omega->vals[idx] + num->rvbOmCnt;

			delta->vals[elemIdx][idx].pib = vXvSparse(lambdaPi, &bOmega);
			if ( num->rvCOmCnt != 0 ) {
				piCrossC = vxMSparse(lambdaPi, &COmega, num->prevCols);
				delta->vals[elemIdx][idx].piC = reduceVector(piCrossC, coord->rvCOmCols, num->rvCOmCnt);
				mem_free(piCrossC);
			}
			else
				delta->vals[elemIdx][idx].piC = NULL;
		}
		mem_free(lambdaPi);
	}

	return 0;
}//END calcDelta()

/* This function allocates memory for an omega structure.  It allocates the memory to structure elements: a dVector to hold an array of
 * observation and the probability associated with it. */
omegaType *newOmega(stocType *stoc, int maxObs, int numStats) {
	omegaType *omega;
	int cnt, i, base, idx;

	omega = (omegaType *) mem_malloc(sizeof(omegaType));
	omega->probs = (dVector) arr_alloc(maxObs, double);
	omega->weights = (iVector) arr_alloc(maxObs, int);
	omega->vals = (dVector *) arr_alloc(maxObs, dVector);
	omega->sampleStats = (dVector *) arr_alloc(numStats, double);
	for ( int n = 0; n < numStats; n++ ) {
		omega->sampleStats[n] = (dVector) arr_alloc(stoc->numOmega+1, double);
	}
	omega->numObs = 0;
	omega->numOmega = stoc->numOmega;
	omega->numStats = numStats;

//	if ( config.ALGO_TYPE == SD ) {
//		return omega;
//	}
//
//	/* Build the entire omega structure if full sample is being used. */
//	if ( strstr(stoc->type, "BLOCKS") != NULL ) {
//		if ( (omega->numObs = stoc->numVals[0]) <= config.MAX_OBS) {
//			omega->vals = (dVector *) mem_realloc(omega->vals, omega->numObs*sizeof(dVector));
//			omega->probs = (dVector) mem_realloc(omega->probs, omega->numObs*sizeof(double));
//			for ( cnt = 0; cnt < omega->numObs; cnt++) {
//				omega->probs[cnt]= stoc->probs[0][cnt];
//				if ( !(omega->vals[cnt] = (dVector) arr_alloc(omega->numOmega+1, double)) )
//					errMsg("allocation", "updateOmega", "omega->vals[cnt]", 0);
//				for (i = 0; i < omega->numOmega; i++)
//					omega->vals[cnt][i+1] = stoc->vals[i][cnt] - stoc->mean[i];
//				omega->vals[cnt][0] = oneNorm(omega->vals[cnt]+1, omega->numOmega);
//			}
//		}
//		else {
//			omega->numObs = config.MAX_OBS;
//		}
//	}
//	else if ( strstr(stoc->type, "INDEP") != NULL ) {
//		omega->numObs = 1; i = 0;
//		while ( i < stoc->numOmega ) {
//			omega->numObs *= stoc->numVals[i];
//			if (omega->numObs > config.MAX_OBS) {
//				omega->numObs = config.MAX_OBS;
//				break;
//			}
//			i++;
//		}
//
//		omega->vals = (dVector *) mem_realloc(omega->vals, omega->numObs*sizeof(dVector));
//		omega->probs = (dVector) mem_realloc(omega->probs, omega->numObs*sizeof(double));
//		for ( cnt = 0; cnt < omega->numObs; cnt++) {
//			if ( !(omega->vals[cnt] = (dVector) arr_alloc(omega->numOmega+1, double)) )
//				errMsg("allocation", "updateOmega", "omega->vals[cnt]", 0);
//			omega->probs[cnt] = 1; base = omega->numObs;
//			for ( i = 0; i < omega->numOmega; i++ ) {
//				base /= stoc->numVals[i];
//				idx = (int)((double) cnt / (double) base) % stoc->numVals[i];
//				omega->vals[cnt][i+1] = stoc->vals[i][idx] - stoc->mean[i];
//				omega->probs[cnt] *= stoc->probs[i][idx];
//			}
//		}
//	}
//	else {
//		omega->numObs = config.MAX_OBS;
//	}

	return omega;
}//END newOmega()

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

void freeLambdaType(lambdaType *lambda, bool partial) {
	int n;

	if (lambda) {
		if (lambda->vals) {
			for ( n = 0; n < lambda->cnt; n++ )
				if (lambda->vals[n]) mem_free(lambda->vals[n]);
			if ( partial ) {
				lambda->cnt = 0;
				return;
			}
			mem_free(lambda->vals);
		}
		mem_free(lambda);
	}

}//END freeLambdaType()

void freeSigmaType(sigmaType *sigma, bool partial) {
	int n;

	if (sigma) {
		for ( n = 0; n < sigma->cnt; n++ )
			if (sigma->vals[n].piC) mem_free(sigma->vals[n].piC);
		if ( partial ) {
			sigma->cnt = 0;
			return;
		}
		if (sigma->lambdaIdx) mem_free(sigma->lambdaIdx);
		if (sigma->vals) mem_free(sigma->vals);
		if (sigma->ck) mem_free(sigma->ck);
		mem_free(sigma);
	}

}//END freeSigmaType()

void freeDeltaType (deltaType *delta, int numDeltaRows, int omegaCnt, bool partial) {
	int n, m;

	if (delta) {
		if (delta->vals) {
			for ( n = 0; n < numDeltaRows; n++ ) {
				if (delta->vals[n]) {
					for ( m = 0; m < omegaCnt; m++ )
						if (delta->vals[n][m].piC)
							mem_free(delta->vals[n][m].piC);
					mem_free(delta->vals[n]);
				}
			}
			if ( partial )
				return;
			mem_free(delta->vals);
		}
		mem_free(delta);
	}

}//END freeDeltaType()

void freeOmegaType(omegaType *omega, bool partial) {

	if ( omega->vals ) {
		for ( int n = 0; n < omega->numObs; n++ )
			if ( omega->vals[n] )
				mem_free(omega->vals[n]);
		if ( partial ) {
			omega->numObs = 0;
			omega->numObs = 0;
			return;
		}
		mem_free(omega->vals);
	}
	if ( omega->probs ) mem_free(omega->probs);
	if ( omega->weights ) mem_free(omega->weights);
	if ( omega->sampleStats) {
		for ( int n = 0; n < omega->numStats; n++ ) {
			if ( omega->sampleStats[n]) mem_free(omega->sampleStats[n]);
		}
		mem_free(omega->sampleStats);
	}

	mem_free(omega);

}//END freeOmegaType()
