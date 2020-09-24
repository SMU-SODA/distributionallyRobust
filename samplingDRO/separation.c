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
extern cString outputDir;

int obtainProbDist(oneProblem *sep, dVector stocMean, omegaType *omega, dVector spObj, int obsStar, bool newOmegaFlag) {
	int status;
	dVector tempProbs;

	/* Update the distribution separation problem */
	if ( updtDistSepProb_MM(sep, stocMean, omega, spObj, obsStar, newOmegaFlag, 1) ) {
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
	/* print the extremal probability distribution to the problem .*/
	FILE *sFile = openFile(outputDir, "extremeProb.csv", "a");
	printVector(omega->probs-1, omega->cnt, sFile);
#endif

	return 0;
}//END obtainProbDist();

int updtDistSepProb_MM(oneProblem *sep, dVector stocMean, omegaType *omega, dVector spObj, int obsStar, bool newOmegaFlag, int numMoments) {
	iVector indices;
	dVector rhsx;
	int idx = 0;

	if ( config.SAMPLING_TYPE == 2 && newOmegaFlag ) {
		/* Add a new column to the distribution separation problem corresponding to the latest observation */
		iVector cmatind;
		dVector cmatval;

		int nzcnt = 2*numMoments*omega->numRV + 1;
		double bdu = 1.0, bdl = 0.0;
		char colname[NAMESIZE];	sprintf(colname,"obs[%d]", obsStar);
		int matbeg = 0;

		cmatind = (iVector) arr_alloc(nzcnt, int);
		cmatval = (dVector) arr_alloc(nzcnt, double);
		for  ( int m = 0; m < numMoments; m++ ) {
			for ( int rv = 0; rv < omega->numRV; rv++ ) {
				cmatind[idx] = idx; idx++;
				cmatval[idx] = stocMean[rv] + omega->vals[obsStar][rv+1];
				cmatind[idx] = idx; idx++;
				cmatval[idx] = -(stocMean[rv] + omega->vals[obsStar][rv+1]);
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
		indices = (iVector) arr_alloc(2*numMoments*omega->numRV, int);
		rhsx    = (dVector) arr_alloc(2*numMoments*omega->numRV+1, double);
		idx = 0;
		for ( int rv = 0; rv < omega->numRV; rv++ ) {
			for  ( int m = 0; m < numMoments; m++ ) {
				indices[idx] = idx;
				rhsx[idx++] = (stocMean[rv] + omega->sampleMean[rv+1])*(1 + config.DRO_PARAM);
				indices[idx] = idx;
				rhsx[idx++] = -(stocMean[rv] + omega->sampleMean[rv+1])*(1 - config.DRO_PARAM);
			}
		}

		if ( changeRHS(sep->lp, idx, indices, rhsx) ) {
			errMsg("solver", "updtDistSepProb", "failed to change the cost coefficients in the solver", 0);
			return 1;
		}

		mem_free(indices);
		mem_free(rhsx);
	}

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

oneProblem *newDistSepProb(dVector stocMean, omegaType *omega) {
	oneProblem *dist = NULL;

	if ( config.DRO_TYPE == RISK_NEUTRAL ) {
		return dist;
	}
	else if ( config.DRO_TYPE == MOMENT_MATCHING ) {
		dist = newDistSepProb_MM(stocMean, omega, 1);
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

oneProblem *newDistSepProb_MM(dVector stocMean, omegaType *omega, int numMoments) {
	oneProblem *dist;
	char tempName[NAMESIZE];
	int offset;

	if (!(dist = (oneProblem *) mem_malloc (sizeof(oneProblem))))
		errMsg("Memory allocation", "newDistSepProb_MM", "failed to allocate memory to distSepProb", 0);

	/* Initialize essential elements */
	dist->type   = PROB_LP;
	dist->objsen = -1;                 				/* sense of the objective: 1 for minimization and -1 for maximization */
	dist->mar 	 = 2*numMoments*omega->numRV + 1;  	/* number of rows is equal to number of moments times number of random variables plus one */
	dist->mac    = omega->cnt;						/* number of columns is equal to the number of observations */
	dist->numInt = 0;                 				/* number of integer variables in the problem  */

	dist->marsz   = dist->mar;
	dist->macsz   = dist->mac;
	dist->matsz   = dist->mar*dist->mac;
	dist->rstorsz = dist->mar*NAMESIZE;
	dist->cstorsz = dist->mac*NAMESIZE;
	dist->numnz   = 0;

	/* Allocate memory to the information whose type is cString */
	if (!(dist->name = (cString) arr_alloc(NAMESIZE, char)))
		errMsg("Allocation", "newDistSepProb_MM", "Fail to allocate memory to dist->name",0);
	if (!(dist->senx = (cString) arr_alloc(dist->marsz,char)))
		errMsg("Allocation", "newDistSepProb_MM", "Fail to allocate memory to dist->senx",0);
	if (!(dist->ctype = (cString) arr_alloc(dist->macsz,char)))
		errMsg("Allocation", "newDistSepProb_MM", "Fail to allocate memory to dist->ctype",0);
	if (!(dist->objname = (cString) arr_alloc(NAMESIZE,char)))
		errMsg("Allocation", "newDistSepProb_MM", "Fail to allocate memory to dist->objname",0);
	if (!(dist->cname = (cString*) arr_alloc(dist->macsz,cString)))
		errMsg("Allocation", "newDistSepProb_MM", "Fail to allocate memory to dist->cname",0);
	if (!(dist->cstore = (cString) arr_alloc(dist->cstorsz, char)))
		errMsg("Allocation", "newDistSepProb_MM", "Fail to allocate memory to dist->cstore",0);
	if (!(dist->rname = (cString *) arr_alloc(dist->marsz,cString)))
		errMsg("Allocation", "newDistSepProb_MM", "Fail to allocate memory to dist->rname",0);
	if (!(dist->rstore = (cString) arr_alloc(dist->rstorsz, char)))
		errMsg("Allocation", "newDistSepProb_MM", "Fail to allocate memory to dist->rstore",0);

	/* Allocate memory to the information whose type is dVector */
	if (!(dist->objx = (dVector) arr_alloc(dist->macsz, double)))
		errMsg("Allocation", "newDistSepProb_MM", "Fail to allocate memory to dist->objx",0);
	if (!(dist->rhsx = (dVector) arr_alloc(dist->marsz, double)))
		errMsg("Allocation", "newDistSepProb_MM", "Fail to allocate memory to dist->rhsx",0);
	if (!(dist->matval = (dVector) arr_alloc(dist->matsz, double)))
		errMsg("allocation", "newDistSepProb_MM", "dist->matval",0);
	if (!(dist->bdl = (dVector) arr_alloc(dist->macsz, double)))
		errMsg("allocation", "newDistSepProb_MM", "dist->bdl",0);
	if (!(dist->bdu = (dVector) arr_alloc(dist->macsz, double)))
		errMsg("allocation", "newDistSepProb_MM", "dist->bdu",0);

	/* Allocate memory to the information whose type is iVector */
	if (!(dist->matbeg = (iVector) arr_alloc(dist->macsz, int)))
		errMsg("allocation", "newDistSepProb_MM", "dist->matbeg",0);
	if (!(dist->matcnt = (iVector) arr_alloc(dist->macsz, int)))
		errMsg("allocation", "newDistSepProb_MM", "dist->matcnt",0);
	if (!(dist->matind = (iVector) arr_alloc(dist->matsz, int)))
		errMsg("allocation", "newDistSepProb_MM", "dist->matind",0);


	strcpy(dist->name, "DistSepProb");           /* Copy problem name */
	strcpy(dist->objname, "DistSep_Obj");     /* Copy objective name */

	/* Add columns to the problem: one for every scenario */
	offset = 0;
	for ( int obs = 0; obs < dist->mac; obs++ ) {
		sprintf(tempName,"obs[%d]", obs);
		strcpy(dist->cstore + offset, tempName);
		dist->cname[obs]  = dist->cstore + offset;
		dist->objx[obs]   = 1.0;
		dist->ctype[obs]  = 'C';
		dist->bdu[obs]    = 1.0;
		dist->bdl[obs]    = 0.0;
		dist->matbeg[obs] = dist->numnz;

		int rowID = 0;
		for  ( int m = 0; m < numMoments; m++ ) {
			for ( int rv = 0; rv < omega->numRV; rv++ ) {
				dist->matind[dist->numnz]   = rowID++;
				dist->matval[dist->numnz++] = stocMean[rv] + omega->vals[obs][rv+1];
				dist->matind[dist->numnz]   = rowID++;
				dist->matval[dist->numnz++] = -(stocMean[rv] + omega->vals[obs][rv+1]);
				dist->matcnt[obs] += 2;
			}
		}

		dist->matind[dist->numnz]   = rowID;
		dist->matval[dist->numnz++] = 1;
		dist->matcnt[obs]++;

		offset += strlen(tempName) + 1;
	}

	/* Add rows to the problem: one for every random variable and moment */
	int j = 0; offset = 0;
	for ( int m = 0; m < numMoments; m++ ) {
		for ( int rv = 0; rv < omega->numRV; rv++ ) {
			sprintf(tempName,"mm_up[%d,%d]", m, rv);
			strcpy(dist->rstore + offset, tempName);
			dist->rname[j] = dist->rstore + offset;
			dist->rhsx[j]  = (stocMean[rv] + omega->sampleMean[rv+1])*(1 + config.DRO_PARAM);
			dist->senx[j++] = 'L';
			offset += strlen(tempName) + 1;

			sprintf(tempName,"mm_dw[%d,%d]", m, rv);
			strcpy(dist->rstore + offset, tempName);
			dist->rname[j]  = dist->rstore + offset;
			dist->rhsx[j]  = -(stocMean[rv] + omega->sampleMean[rv+1])*(1 - config.DRO_PARAM);
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

#if defined(SETUP_CHECK)
	if ( writeProblem(dist->lp, "newDistSep.lp") ) {
		errMsg("solver", "newDistSepProb_MM", "failed to write distribution separation problem to file", 0);
		return NULL;
	}
#endif

	return dist;
}//END newDistProb()
