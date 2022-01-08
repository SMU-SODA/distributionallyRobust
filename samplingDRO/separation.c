/*
 * separation.c
 *
 *  Created on: Jul 22, 2019
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "samplingDRO.h"

extern configType config;

#if defined(SEP_CHECK)
extern cString outputDir;
#endif

int obtainProbDist(oneProblem *sep, dVector stocMean, omegaType *omega, dVector spObj, int obsStar, bool newOmegaFlag) {
	int status;
	dVector tempProbs;

	/* Update the distribution separation problem */
	if ( updtDistSepProb_MM(sep, stocMean, omega, spObj, obsStar, newOmegaFlag) ) {
		errMsg("separation", "updtDistSepProb_MM", "failed to update the distribution separation problem", 0);
		return 1;
	}

	/* Solve the distribution separation problem */
	if ( solveProblem(sep->lp, sep->name, sep->type, &status) ) {
		if ( status == STAT_INFEASIBLE ) {
			printf("Distribution separation problem is infeasible.\n");
			return 1;
		}
		else {
			errMsg("algorithm", "obtainProbDist", "failed to solve subproblem in solver", 0);
			return 1;
		}
	}

	/* Obtain the extremal probability distribution */
	tempProbs = (dVector) arr_alloc(omega->cnt + 1, double);
	if ( getPrimal(sep->lp, tempProbs, omega->cnt) ) {
		errMsg("algorithm", "obtainProbDist", "failed to obtain the primal solution for master", 0);
		return 1;
	}
	for ( int obs = 0; obs < omega->cnt; obs++ ) {
		omega->probs[obs] = tempProbs[obs+1];
	}
	mem_free(tempProbs);

#if defined(SEP_CHECK)
	/* print the extremal probability distribution to the problem.*/
	FILE *sFile = openFile(outputDir, "extremeProb.csv", "a");
	printVector(omega->probs-1, omega->cnt, sFile);
	fclose(sFile);
#endif

	return 0;
}//END obtainProbDist();

int updtDistSepProb_MM(oneProblem *sep, dVector stocMean, omegaType *omega, dVector spObj, int obsStar, bool newOmegaFlag) {
	iVector indices;
	dVector rhsx;
	int idx, numMoments = config.DRO_PARAM_1;

	if ( config.ALGO_TYPE == SD && newOmegaFlag ) {
		/* Add a new column to the distribution separation problem corresponding to the latest observation */
		iVector cmatind;
		dVector cmatval;

		int nzcnt = 2*numMoments*omega->numOmega + 1;
		double bdu = 1.0, bdl = 0.0;
		char colname[NAMESIZE];	sprintf(colname,"obs[%d]", obsStar);
		int matbeg = 0;

		cmatind = (iVector) arr_alloc(nzcnt, int);
		cmatval = (dVector) arr_alloc(nzcnt, double);
		idx = 0;
		for ( int rv = 1; rv <= omega->numOmega; rv++ ) {
			for  ( int m = 0; m < numMoments; m++ ) {
				cmatind[idx] = idx;
				cmatval[idx++] = pow(stocMean[rv] + omega->vals[obsStar][rv], m+1);
				cmatind[idx] = idx;
				cmatval[idx++] = -pow(stocMean[rv] + omega->vals[obsStar][rv], m+1);
			}
		}
		cmatind[idx] = idx;
		cmatval[idx] = 1.0;

		if ( addCol(sep->lp, nzcnt, 1.0, matbeg, cmatind, cmatval, bdu, bdl, colname) ) {
			errMsg("solver", "updtDistSepProb", "failed to add a new column for latest observation", 0);
			return 1;
		}

		mem_free(cmatind); mem_free(cmatval);

		/* Update the right-hand side with sample mean information */
		indices = (iVector) arr_alloc(2*numMoments*omega->numOmega, int);
		rhsx    = (dVector) arr_alloc(2*numMoments*omega->numOmega+1, double);
		idx = 0;
		for ( int rv = 1; rv <= omega->numOmega; rv++ ) {
			for  ( int m = 0; m < numMoments; m++ ) {
				indices[idx] = idx;
				rhsx[idx++] = omega->sampleStats[m][rv]*(1 + config.DRO_PARAM_2);
				indices[idx] = idx;
				rhsx[idx++] = -omega->sampleStats[m][rv]*(1 - config.DRO_PARAM_2);
			}
		}

		if ( changeRHS(sep->lp, idx, indices, rhsx) ) {
			errMsg("solver", "updtDistSepProb", "failed to change the cost coefficients in the solver", 0);
			return 1;
		}

		mem_free(indices);
		mem_free(rhsx);
	}

	writeProblem(sep->lp, "cellSepProb.lp");

	/* Update the objective coefficients with the recourse function estimates */
	indices = (iVector) arr_alloc(omega->cnt, int);
	for ( int n = 0; n < omega->cnt; n++ ) {
		indices[n] = n;
	}

	if ( changeObjx(sep->lp, omega->cnt, indices, spObj) ) {
		errMsg("solver", "updtDistSepProb", "failed to change the cost coefficients in the solver", 0);
		return 1;
	}

#if defined( SEP_CHECK )
	writeProblem(sep->lp, "cellSepProb.lp");
#endif

	mem_free(indices);
	return 0;
}//END updtDistSepProb()

/* Setup a new distribution separation problem. Currently there are two types of ambiguity sets supported by this implementation.
 *
 * 		1. Moment-based ambiguity set
 * 		2. Wasserstein distance-based ambiguity set
 *
 * The distribution separation problem uses the oneProblem structure, but the elements of the structure are not updated when the
 * distribution separation problem gets updated during the course of the algorithm and/or replications. In this sense, only the
 * LP pointer is useful for most purposes.
 */
oneProblem *newDistSepProb(dVector stocMean, omegaType *omega) {
	oneProblem *dist = NULL;

	if ( config.DRO_TYPE == RISK_NEUTRAL ) {
		return dist;
	}
	else if ( config.DRO_TYPE == MOMENT_MATCHING ) {
		/* Obtain the necessary statistics for creating a distribution separation problem */
		refineOmega(omega, stocMean);

		dist = newDistSepProb_MM(stocMean, omega);
		if ( dist == NULL ) {
			errMsg("solve", "newDistSepProb", "unknown distribution type requested in DRO_TYPE", 0);
			return NULL;
		}
	}
	else if ( config.DRO_TYPE == WASSERSTEIN ) {
		dist = newDistSepProb_W(omega);
		if ( dist == NULL ) {
			errMsg("solve", "newDistSepProb", "unknown distribution type requested in DRO_TYPE", 0);
			return NULL;
		}
	}
	else {
		errMsg("solve", "newDistSepProb", "unknown distribution type requested in DRO_TYPE", 0);
		return NULL;
	}

	return dist;
}//END newDistSeparation()

/* This subroutine is used to setup a distribution separation problem with moment-based ambiguity set. */
oneProblem *newDistSepProb_MM(dVector stocMean, omegaType *omega) {
	oneProblem *dist;
	char tempName[NAMESIZE];
	int offset, numMoments = config.DRO_PARAM_1;

	if (!(dist = (oneProblem *) mem_malloc (sizeof(oneProblem))))
		errMsg("Memory allocation", "newDistSepProb_MM", "failed to allocate memory to distSepProb", 0);

	/* Initialize essential elements */
	dist->type   = PROB_LP;
	dist->objsen = -1;                 				/* sense of the objective: 1 for minimization and -1 for maximization */
	dist->mar 	 = 2*numMoments*omega->numOmega + 1;  	/* number of rows is equal to number of moments times number of random variables plus one */
	dist->mac    = omega->cnt;						/* number of columns is equal to the number of observations */
	dist->numInt = 0;                 				/* number of integer variables in the problem  */

	dist->marsz   = dist->mar;
	dist->macsz   = dist->mac;
	dist->matsz   = dist->mar*dist->mac;
	dist->rstorsz = dist->mar*NAMESIZE;
	dist->cstorsz = dist->mac*NAMESIZE;
	dist->numnz   = 0;

	/* Allocate memory to the information whose type is cString */
	dist->name = (cString) arr_alloc(NAMESIZE, char);
	dist->senx = (cString) arr_alloc(dist->marsz,char);
	dist->ctype = (cString) arr_alloc(dist->macsz,char);
	dist->objname = (cString) arr_alloc(NAMESIZE,char);
	dist->cname = (cString*) arr_alloc(dist->macsz,cString);
	dist->cstore = (cString) arr_alloc(dist->cstorsz, char);
	dist->rname = (cString *) arr_alloc(dist->marsz,cString);
	dist->rstore = (cString) arr_alloc(dist->rstorsz, char);

	/* Allocate memory to the information whose type is dVector */
	dist->objx = (dVector) arr_alloc(dist->macsz, double);
	dist->rhsx = (dVector) arr_alloc(dist->marsz, double);
	dist->matval = (dVector) arr_alloc(dist->matsz, double);
	dist->bdl = (dVector) arr_alloc(dist->macsz, double);
	dist->bdu = (dVector) arr_alloc(dist->macsz, double);

	/* Allocate memory to the information whose type is iVector */
	dist->matbeg = (iVector) arr_alloc(dist->macsz, int);
	dist->matcnt = (iVector) arr_alloc(dist->macsz, int);
	dist->matind = (iVector) arr_alloc(dist->matsz, int);

	strcpy(dist->name, "DistSepProb_MM");          	/* Copy problem name */
	strcpy(dist->objname, "DistSep_Obj");     		/* Copy objective name */

	/* Add columns to the problem: one for every scenario */
	offset = 0;
	for ( int obs = 0; obs < dist->mac; obs++ ) {
		sprintf(tempName,"obs[%d]", obs);
		strcpy(dist->cstore + offset, tempName);
		dist->cname[obs]  = dist->cstore + offset;
		offset += strlen(tempName) + 1;
		dist->objx[obs]   = 1.0;
		dist->ctype[obs]  = 'C';
		dist->bdu[obs]    = 1.0;
		dist->bdl[obs]    = 0.0;
		dist->matbeg[obs] = dist->numnz;

		int rowID = 0;
		for ( int rv = 1; rv <= omega->numOmega; rv++ ) {
			for  ( int m = 0; m < numMoments; m++ ) {
				dist->matind[dist->numnz]   = rowID++;
				dist->matval[dist->numnz++] = pow(stocMean[rv] + omega->vals[obs][rv], m+1);
				dist->matind[dist->numnz]   = rowID++;
				dist->matval[dist->numnz++] = -pow(stocMean[rv] + omega->vals[obs][rv], m+1);
				dist->matcnt[obs] += 2;
			}
		}

		dist->matind[dist->numnz]   = rowID;
		dist->matval[dist->numnz++] = 1;
		dist->matcnt[obs]++;
	}

	/* Add rows to the problem: one for every random variable and moment */
	int j = 0; offset = 0;
	for ( int rv = 1; rv <= omega->numOmega; rv++ ) {
		for ( int m = 0; m < numMoments; m++ ) {
			sprintf(tempName,"mm_up[%d,%d]", m, rv-1);
			strcpy(dist->rstore + offset, tempName);
			dist->rname[j] = dist->rstore + offset;
			dist->rhsx[j]  = omega->sampleStats[m][rv]*(1 + config.DRO_PARAM_2);
			dist->senx[j++] = 'L';
			offset += strlen(tempName) + 1;

			sprintf(tempName,"mm_dw[%d,%d]", m, rv-1);
			strcpy(dist->rstore + offset, tempName);
			dist->rname[j]  = dist->rstore + offset;
			dist->rhsx[j]  = -omega->sampleStats[m][rv]*(1 - config.DRO_PARAM_2);
			dist->senx[j++] = 'L';
			offset += strlen(tempName) + 1;
		}
	}
	sprintf(tempName, "prob");
	strcpy(dist->rstore + offset, tempName);
	dist->rname[j] = dist->rstore + offset;
	dist->rhsx[j]  = 1.0;
	dist->senx[j]  = 'E';

	/* Load the copy into CPLEX */
	dist->lp = setupProblem(dist->name, dist->type, dist->mac, dist->mar, dist->objsen, dist->objx, dist->rhsx, dist->senx,
			dist->matbeg, dist->matcnt,dist->matind, dist->matval, dist->bdl, dist->bdu, NULL, dist->cname, dist->rname,
			dist->ctype);
	if ( dist->lp == NULL ) {
		errMsg("Problem Setup", "newDistSepProb_MM", "failed to setup master problem in the solver",0);
		return NULL;
	}

#if defined(SEP_CHECK)
	if ( writeProblem(dist->lp, "newDistSep.lp") ) {
		errMsg("solver", "newDistSepProb_MM", "failed to write distribution separation problem to file", 0);
		return NULL;
	}
#endif

	return dist;
}//END newDistProb()

oneProblem *newDistSepProb_W(omegaType *omega) {
	oneProblem *dist = NULL;
	char tempName[NAMESIZE];
	double **distMatrix;

	/* Calculate the distance matrix between all the observation */
	distMatrix = (double **) arr_alloc(omega->cnt, double *);
	for ( int i = 0; i < omega->cnt; i++ ) {
		distMatrix[i] = (double *) arr_alloc(omega->cnt, double);
		for (int j = i+1; j < omega->cnt; j++ ) {
			distMatrix[i][j] = pNorm(omega->vals[i], omega->vals[j], omega->numOmega, config.DRO_PARAM_1);
		}
	}

	/* Allocate memory */
	dist = (oneProblem *) mem_malloc (sizeof(oneProblem));

	/* Initialize essential elements */
	dist->type   = PROB_LP;
	dist->objsen = -1;                 			/* sense of the objective: 1 for minimization and -1 for maximization */
	dist->mar 	 = 2*omega->cnt + 2;			/* number of rows is equal to number of moments times number of random variables plus one */
	dist->mac    = omega->cnt*(omega->cnt + 1);	/* number of columns is equal to the number of observations */
	dist->numInt = 0;                 			/* number of integer variables in the problem  */

	dist->marsz   = dist->mar;
	dist->macsz   = dist->mac;
	dist->matsz   = dist->mar*dist->mac;
	dist->rstorsz = dist->mar*NAMESIZE;
	dist->cstorsz = dist->mac*NAMESIZE;
	dist->numnz   = 0;

	/* Allocate memory to the information whose type is cString */
	dist->name = (cString) arr_alloc(NAMESIZE, char);
	dist->senx = (cString) arr_alloc(dist->marsz,char);
	dist->ctype = (cString) arr_alloc(dist->macsz,char);
	dist->objname = (cString) arr_alloc(NAMESIZE,char);
	dist->cname = (cString*) arr_alloc(dist->macsz,cString);
	dist->cstore = (cString) arr_alloc(dist->cstorsz, char);
	dist->rname = (cString *) arr_alloc(dist->marsz,cString);
	dist->rstore = (cString) arr_alloc(dist->rstorsz, char);

	/* Allocate memory to the information whose type is dVector */
	dist->objx = (dVector) arr_alloc(dist->macsz, double);
	dist->rhsx = (dVector) arr_alloc(dist->marsz, double);
	dist->matval = (dVector) arr_alloc(dist->matsz, double);
	dist->bdl = (dVector) arr_alloc(dist->macsz, double);
	dist->bdu = (dVector) arr_alloc(dist->macsz, double);

	/* Allocate memory to the information whose type is iVector */
	dist->matbeg = (iVector) arr_alloc(dist->macsz, int);
	dist->matcnt = (iVector) arr_alloc(dist->macsz, int);
	dist->matind = (iVector) arr_alloc(dist->matsz, int);

	strcpy(dist->name, "DistSepProb_W");          	/* Copy problem name */
	strcpy(dist->objname, "DistSep_Obj");     		/* Copy objective name */

	/* Add columns to the problem: one for every scenario */
	/* a. marginal variables (p) */
	int offset = 0, rowID = 0, colID = 0;
	for ( int obs = 0; obs < omega->cnt; obs++ ) {

		sprintf(tempName,"prob[%d]", obs);
		strcpy(dist->cstore + offset, tempName);
		dist->cname[colID]  = dist->cstore + offset;
		offset += strlen(tempName) + 1;
		dist->objx[colID]   = 1.0;
		dist->ctype[colID]  = 'C';
		dist->bdu[colID]    = 1.0;
		dist->bdl[colID]    = 0.0;
		dist->matbeg[colID] = dist->numnz;

		dist->matind[dist->numnz] = 0;
		dist->matval[dist->numnz++] = 1.0;
		dist->matind[dist->numnz] = obs+1;
		dist->matval[dist->numnz++] = -1.0;

		dist->matcnt[colID++] += 2;
	}

	/* b. joint probability variables (eta) */
	for ( int i = 0; i < omega->cnt; i++ ) {
		for ( int j = 0; j < omega->cnt; j++ ) {
			if ( i != j ) {
				sprintf(tempName,"eta[%d][%d]", i, j);
				strcpy(dist->cstore + offset, tempName);
				dist->cname[colID]  = dist->cstore + offset;
				offset += strlen(tempName) + 1;
				dist->objx[colID]   = 0.0;
				dist->ctype[colID]  = 'C';
				dist->bdu[colID]    = INFINITY;
				dist->bdl[colID]    = 0.0;
				dist->matbeg[colID] = dist->numnz;

				dist->matind[dist->numnz]   = i+1;
				dist->matval[dist->numnz++] = 1.0;

				dist->matind[dist->numnz]   = omega->cnt+j+1;
				dist->matval[dist->numnz++] = 1.0;

				dist->matind[dist->numnz]   = 2*omega->cnt+1;
				dist->matval[dist->numnz++] = distMatrix[i][j];

				dist->matcnt[colID++] += 3;
			}
		}
	}

	/* Add rows to the problem: one for every random variable and moment */
	/* a. probability sums to one */
	rowID = 0; offset = 0;
	sprintf(tempName,"probSum");
	strcpy(dist->rstore + offset, tempName);
	dist->rname[rowID] = dist->rstore + offset;
	dist->rhsx[rowID]  = 1.0;
	dist->senx[rowID++] = 'E';
	offset += strlen(tempName) + 1;

	/* b. row sum for joint probability */
	for ( int rv = 1; rv <= omega->cnt; rv++ ) {
		sprintf(tempName,"rowSum[%d]", rv-1);
		strcpy(dist->rstore + offset, tempName);
		dist->rname[rowID] = dist->rstore + offset;
		dist->rhsx[rowID]  = 0.0;
		dist->senx[rowID++] = 'E';
		offset += strlen(tempName) + 1;
	}

	/* c. column sum for joint probability */
	for ( int rv = 1; rv <= omega->cnt; rv++ ) {
		sprintf(tempName,"colSum[%d]", rv-1);
		strcpy(dist->rstore + offset, tempName);
		dist->rname[rowID] = dist->rstore + offset;
		dist->rhsx[rowID]  = (double) omega->weights[rv-1]/(double) omega->numObs;
		dist->senx[rowID++] = 'E';
		offset += strlen(tempName) + 1;
	}

	/* d. limit on the distance */
	sprintf(tempName,"distLim");
	strcpy(dist->rstore + offset, tempName);
	dist->rname[rowID] = dist->rstore + offset;
	dist->rhsx[rowID]  = config.DRO_PARAM_2;
	dist->senx[rowID++] = 'L';

	/* Load the copy into CPLEX dist->mac, dist->mar,*/
	dist->lp = setupProblem(dist->name, dist->type, colID, rowID, dist->objsen, dist->objx, dist->rhsx, dist->senx,
			dist->matbeg, dist->matcnt,dist->matind, dist->matval, dist->bdl, dist->bdu, NULL, dist->cname, dist->rname,
			dist->ctype);
	if ( dist->lp == NULL ) {
		errMsg("Problem Setup", "newDistSepProb_MM", "failed to setup master problem in the solver",0);
		return NULL;
	}

#if defined(SEP_CHECK)
	if ( writeProblem(dist->lp, "newDistSep.lp") ) {
		errMsg("solver", "newDistSepProb_MM", "failed to write distribution separation problem to file", 0);
		return NULL;
	}
#endif

	if ( distMatrix ) {
		for ( int i = 0; i < omega->cnt; i++ ) {
			mem_free(distMatrix[i]);
		}
		mem_free(distMatrix);
	}

	return dist;
}//END newDistSepProb_W()

int cleanDistSepProb(oneProblem *sep, dVector stocMean, omegaType *omega, int numMoments) {
	dVector rhsx;
	iVector indices;
	int idx;

	/* Remove all the columns in the separation problem. For sequential sampling, a clean distribution separation problem is
	 * needed. For external sampling, the constraint coefficients in all columns need to be updated. Therefore, removing all
	 * columns in efficient.
	 */
	{
		int mac = getNumCols(sep->lp);
		if (  removeColumns(sep->lp, 0, mac-1) ) {
			errMsg("solver", "cleanDistSepProb", "failed to remove a column from the distribution separation problem problem", 0);
			return 1;
		}
	}

	if ( config.ALGO_TYPE != SD ) {
		/* Add columns corresponding to observations in the current replications. */
		iVector cmatind;
		dVector cmatval;
		int nzcnt = 2*numMoments*omega->numOmega + 1;
		double bdu = 1.0, bdl = 0.0;
		char colname[NAMESIZE];
		int matbeg = 0;

		cmatind = (iVector) arr_alloc(nzcnt, int);
		cmatval = (dVector) arr_alloc(nzcnt, double);

		for ( int obs = 0; obs < omega->cnt; obs++ ) {
			sprintf(colname,"obs[%d]", obs);

			idx = 0;
			for ( int rv = 1; rv <= omega->numOmega; rv++ ) {
				for  ( int m = 0; m < numMoments; m++ ) {
					cmatind[idx] = idx;
					cmatval[idx++] = pow(stocMean[rv] + omega->vals[obs][rv], m+1);
					cmatind[idx] = idx;
					cmatval[idx++] = -pow(stocMean[rv] + omega->vals[obs][rv], m+1);
				}
			}
			cmatind[idx] = idx;
			cmatval[idx] = 1.0;

			if ( addCol(sep->lp, nzcnt, 1.0, matbeg, cmatind, cmatval, bdu, bdl, colname) ) {
				errMsg("solver", "updtDistSepProb", "failed to add a new column for latest observation", 0);
				return 1;
			}
		}

		mem_free(cmatind); mem_free(cmatval);

		/* The only update necessary is in the right hand side of the constraints */
		indices = (iVector) arr_alloc(2*numMoments*omega->numOmega, int);
		rhsx    = (dVector) arr_alloc(2*numMoments*omega->numOmega+1, double);
		idx = 0;
		for ( int rv = 1; rv <= omega->numOmega; rv++ ) {
			for  ( int m = 0; m < numMoments; m++ ) {
				indices[idx] = idx;
				rhsx[idx++] = omega->sampleStats[m][rv]*(1 + config.DRO_PARAM_2);
				indices[idx] = idx;
				rhsx[idx++] = -omega->sampleStats[m][rv]*(1 - config.DRO_PARAM_2);
			}
		}

		if ( changeRHS(sep->lp, idx, indices, rhsx) ) {
			errMsg("solver", "cleanDistSepProb", "failed to change the cost coefficients in the solver", 0);
			return 1;
		}

		mem_free(indices);
		mem_free(rhsx);
	}

#if defined(SEP_CHECK)
	if ( writeProblem(sep->lp, "cleanDistSep.lp") ) {
		errMsg("solver", "cleanDistSepProb", "failed to write distribution separation problem to file", 0);
		return 1;
	}
#endif

	return 0;
}//END cleanDistSepProb()
