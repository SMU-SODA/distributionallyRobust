/*
 * reform.c
 *
 *  Created on: Sep 29, 2021
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "samplingDRO.h"
configType config;

int solveReformulation(stocType *stoc, probType **prob, cellType **cell);
oneProblem *setupMomentMatchReform (probType **prob, omegaType *omega);
oneProblem *setupWassersteinInfReform (probType **prob, omegaType *omega);

int solveReformulation(stocType *stoc, probType **prob, cellType **cell) {
	int status;

	(*cell) = createEmptyCell();

	/* Simulate the observations */
	(*cell)->omega  = newOmega(stoc, config.MAX_OBS, 0);
	(*cell)->omega->numObs = config.MAX_OBS;
	setupSAA(stoc, NULL, &config.RUN_SEED[0], &(*cell)->omega->vals, (*cell)->omega->probs, (*cell)->omega->weights,
			&(*cell)->omega->numObs, config.TOLERANCE);

	/* Setup the reformulated problem */
	if ( config.DRO_TYPE == RISK_NEUTRAL ) {
		errMsg("solve", "solveReformulation", "no subroutine to setup extensive scenario formulation", 0);
		return 1;
	}
	else if ( config.DRO_TYPE == MOMENT_MATCHING ) {
		(*cell)->master = setupMomentMatchReform(prob, (*cell)->omega);
		if ( (*cell)->master == NULL ) {
			errMsg("solve", "solveReformulation", "unable to setup moment matching-based reformulation", 0);
			return 1;
		}
	}
	else if ( config.DRO_TYPE == WASSERSTEIN ) {

		if ( config.DRO_PARAM_1 == 1) {
			(*cell)->master = setupWassersteinOneReform (prob, (*cell)->omega);
			if ( (*cell)->master == NULL ) {
				errMsg("solve", "solveReformulation", "unable to setup 1-Wasserstein reformulation", 0);
				return 1;
			}
		}
		else if ( config.DRO_PARAM_1 < 0 ) {
			(*cell)->master = setupWassersteinInfReform (prob, (*cell)->omega);
			if ( (*cell)->master == NULL ) {
				errMsg("solve", "solveReformulation", "unable to setup Inf-Wasserstein reformulation", 0);
				return 1;
			}
		}
		else {
			errMsg("solve", "solveReformulation", "reformulation for the chosen Wasserstein metric is not supported", 0);
			return 1;
		}
	}
	else {
		errMsg("solve", "solveReformulation", "unknown distribution type requested in DRO_TYPE", 0);
		return 1;
	}

	/* Solve the reformulated problem */
	if ( solveProblem((*cell)->master->lp, (*cell)->master->name, config.MASTER_TYPE, &status) ) {
		writeProblem((*cell)->master->lp, "error.lp");
		errMsg("algorithm", "solveReformulation", "failed to solve the master problem", 0);
		return 1;
	}

	/* Obtain the primal solution and the objective function value. */
	(*cell)->candidX = (dVector) arr_alloc(prob[0]->num->cols+1, double);
	if ( getPrimal((*cell)->master->lp, (*cell)->candidX, prob[0]->num->cols) ) {
		errMsg("algorithm", "solveReformulation", "failed to obtain the primal solution for master", 0);
		return 1;
	}
	(*cell)->incumbEst = getObjective((*cell)->master->lp, PROB_LP);

	for ( int i = 0; i < (*cell)->master->mac; i++ )
		mem_free((*cell)->master->cname[i]);
	for ( int i = 0; i < (*cell)->master->mar; i++ )
		mem_free((*cell)->master->rname[i]);
	return 0;
}//END solveReformulation()

oneProblem *setupMomentMatchReform (probType **prob, omegaType *omega) {
	oneProblem *reform = NULL;
	char tempName[NAMESIZE];

	reform = (oneProblem *) mem_malloc (sizeof(oneProblem));
	reform->type = PROB_LP;
	reform->objsen  = 1;

	reform->mar 	= prob[0]->sp->mar + omega->numObs*prob[1]->sp->mar + omega->numObs;				/* number of rows */
	reform->mac 	= prob[0]->sp->mac + omega->numObs*prob[1]->sp->mac + 2*omega->numObs + 1;			/* number of columns + etas */
	reform->numInt 	= prob[0]->sp->numInt + omega->numObs*prob[1]->sp->numInt;    						/* number of integer variables in the problem  */
	reform->numnz 	= prob[0]->sp->numnz + omega->numObs*(prob[1]->sp->numnz + prob[1]->Cbar->cnt)
											+ 3*omega->numObs;											/* number of non-zero elements in constraint matrix */
	reform->matsz 	= prob[0]->sp->matsz + omega->numObs*(prob[1]->sp->matsz + prob[1]->Cbar->cnt)
					+ omega->numObs*(prob[0]->sp->mac + prob[1]->sp->mac + 3);							/* extended matrix size */
	reform->marsz 	= prob[0]->sp->marsz + omega->numObs*prob[1]->sp->marsz + omega->numObs;      		/* extended row size */
	reform->macsz 	= prob[0]->sp->macsz + omega->numObs*prob[1]->sp->macsz + 2*omega->numObs + 1; 		/* extended column size */

	reform->rstorsz = reform->cstorsz = 0;

	/* Allocate memory to the information whose type is cString */
	reform->name = (cString) arr_alloc(NAMESIZE, char);
	reform->senx = (cString) arr_alloc(reform->marsz,char);
	reform->ctype = (cString) arr_alloc(reform->macsz,char);
	reform->objname = (cString) arr_alloc(NAMESIZE,char);
	reform->cname = (cString*) arr_alloc(reform->macsz, cString);
	reform->rname = (cString *) arr_alloc(reform->marsz, cString);
	reform->cstore = NULL;
	reform->rstore = NULL;

	/* Allocate memory to the information whose type is dVector */
	reform->objx = (dVector) arr_alloc(reform->macsz, double);
	reform->rhsx = (dVector) arr_alloc(reform->marsz, double);
	reform->matval = (dVector) arr_alloc(reform->matsz, double);
	reform->bdl = (dVector) arr_alloc(reform->macsz, double);
	reform->bdu = (dVector) arr_alloc(reform->macsz, double);

	/* Allocate memory to the information whose type is iVector */
	reform->matbeg = (iVector) arr_alloc(reform->macsz, int);
	reform->matcnt = (iVector) arr_alloc(reform->macsz, int);
	reform->matind = (iVector) arr_alloc(reform->matsz, int);

	strcpy(reform->name, "momentMatchReform");	/* Copy problem name */
	strcpy(reform->objname, "MomentObj");      	/* Copy objective name */

	/* First-stage problem */
	int cnt = 0;
	/* Copy all the column information from the original first-stage problem */
	for (int j = 0; j < prob[0]->sp->mac; j++) {
		reform->objx[j]   = prob[0]->sp->objx[j];					/* Copy objective function coefficients */
		reform->ctype[j]  = prob[0]->sp->ctype[j];					/* Copy the decision variable type */
		reform->bdu[j] 	  = prob[0]->sp->bdu[j];					/* Copy the upper bound and lower bound */
		reform->bdl[j] 	  = prob[0]->sp->bdl[j];
		reform->matbeg[j] = cnt;									/* Copy the master sparse matrix beginning position of each column */
		reform->matcnt[j] = 0;
		for (int idx = prob[0]->sp->matbeg[j]; idx < prob[0]->sp->matbeg[j] + prob[0]->sp->matcnt[j]; idx++) {
			reform->matval[cnt] = prob[0]->sp->matval[idx];			/* Copy the non-zero coefficient */
			reform->matind[cnt] = prob[0]->sp->matind[idx];			/* Copy the row entry of the non-zero elements */
			reform->matcnt[j]++;
			cnt++;
		}

		/* Loop through the transfer matrix to introduce the first-stage variables for second-stage constraints for all observations */
		for (int idx = 0; idx < prob[1]->Cbar->cnt; idx++ ) {
			if ( prob[1]->Cbar->col[idx+1] == j + 1 ) {
				for ( int obs = 0; obs < omega->numObs; obs++ ) {
					reform->matval[cnt] = prob[1]->Cbar->val[idx+1];
					reform->matind[cnt] = prob[0]->sp->mar + obs*prob[1]->sp->mar + prob[1]->Cbar->row[idx+1]-1;
					reform->matcnt[j]++;
					cnt++;
				}
			}
		}

		/* Add the column coefficients for the additional constraints that ensure deviation. */
		for ( int obs = 0; obs < omega->numObs; obs++ ) {
			reform->matval[cnt] = prob[0]->sp->objx[j];
			reform->matind[cnt] = prob[0]->sp->mar + omega->numObs*prob[1]->sp->mar + obs;
			reform->matcnt[j]++;
			cnt++;
		}

		reform->cname[j] = (cString) arr_alloc(NAMESIZE, char);
		strcpy(reform->cname[j], prob[0]->sp->cname[j]);
	}

	/* Copy all information concerning rows of first-stage */
	for (int r = 0; r < prob[0]->sp->mar; r++) {
		reform->rhsx[r]  = prob[0]->sp->rhsx[r]; 				/* Copy the right hand side value */
		reform->senx[r]  = prob[0]->sp->senx[r];				/* Copy the constraint sense */

		reform->rname[r] = (cString) arr_alloc(NAMESIZE, char);
		strcpy(reform->rname[r], prob[0]->sp->rname[r]);
	}

	/* Second-stage problem */
	/* Loop through all the observations to add the variables and constraints for each one. */
	for ( int obs = 0; obs < omega->numObs; obs++) {
		/* Copy all the column information from the original second-stage problem */
		int k = prob[0]->sp->mac + obs*prob[1]->sp->mac;
		for (int j = 0; j < prob[1]->sp->mac; j++) {
			reform->objx[k+j] 	= prob[1]->sp->objx[j]/(double) omega->numObs;					/* Copy objective function coefficients */
			reform->ctype[k+j] 	= prob[1]->sp->ctype[j];				/* Copy the decision variable type */
			reform->bdu[k+j] 	= prob[1]->sp->bdu[j];					/* Copy the upper bound and lower bound */
			reform->bdl[k+j] 	= prob[1]->sp->bdl[j];
			reform->matbeg[k+j] = cnt;									/* Copy the master sparse matrix beginning position of each column */
			reform->matcnt[k+j] = 0;
			for (int idx = prob[1]->sp->matbeg[j]; idx < prob[1]->sp->matbeg[j] + prob[1]->sp->matcnt[j]; idx++) {
				reform->matval[cnt] = prob[1]->sp->matval[idx];		/* Copy the non-zero coefficient */
				reform->matind[cnt] = prob[0]->sp->mar + obs*prob[1]->sp->mar + prob[1]->sp->matind[idx];		/* Copy the row entry of the non-zero elements */
				reform->matcnt[k+j]++;
				cnt++;
			}

			/* Add the coefficient of the column to the additional moment constraint */
			reform->matval[cnt] = prob[1]->sp->objx[j];
			reform->matind[cnt] = prob[0]->sp->mar + omega->numObs*prob[1]->sp->mar + obs;
			reform->matcnt[k+j]++;
			cnt++;

			reform->cname[k+j] = (cString) arr_alloc(NAMESIZE, char);
			sprintf(tempName, "%s_%d", prob[1]->sp->cname[j], obs);
			strcpy(reform->cname[k+j], tempName);
		}

		/* Copy all information concerning rows of second-stage */
		k = prob[0]->sp->mar + obs*prob[1]->sp->mar;
		for (int r = 0; r < prob[1]->sp->mar; r++) {
			reform->rhsx[k+r]  = prob[1]->sp->rhsx[r]; 				/* Copy the right hand side value */
			reform->senx[k+r]  = prob[1]->sp->senx[r];				/* Copy the constraint sense */

			reform->rname[k+r] = (cString) arr_alloc(NAMESIZE, char);
			sprintf(tempName, "%s_%d", prob[1]->sp->rname[r], obs);
			strcpy(reform->rname[k+r], tempName);
		}
	}

	/* Add the additional lambda variable */
	int k = prob[0]->sp->mac + omega->numObs*prob[1]->sp->mac;
	reform->objx[k]  = 0.0;
	reform->ctype[k] = 'C';
	reform->bdu[k] 	 = INFINITY;
	reform->bdl[k] 	 = 0.0;
	reform->matbeg[k] = cnt;
	reform->matcnt[k] = 0;
	reform->cname[k] = (cString) arr_alloc(NAMESIZE, char);	strcpy(reform->cname[k], "lambda");
	for ( int obs = 0; obs < omega->numObs; obs++) {
		reform->matval[cnt] = -8000;
		reform->matind[cnt] = prob[0]->sp->mar + omega->numObs*prob[1]->sp->mar + obs;
		reform->matcnt[k]++;
		cnt++;
	}

	k++;
	for ( int obs = 0; obs < omega->numObs; obs++) {
		/* Add the additional uPlus and uMinus variables */
		reform->objx[k]  = config.DRO_PARAM_2;  		reform->objx[k+1] = config.DRO_PARAM_2;
		reform->ctype[k] = 'C';			reform->ctype[k+1]  = 'C';
		reform->bdu[k] 	 = INFINITY;	reform->bdu[k+1]    = INFINITY;
		reform->bdl[k] 	 = 0.0;			reform->bdl[k+1]    = 0.0;
		reform->matbeg[k] = cnt;		reform->matbeg[k+1] = cnt+1;
		reform->matcnt[k] = 1;			reform->matcnt[k+1] = 1;

		reform->cname[k] = (cString) arr_alloc(NAMESIZE, char);
		sprintf(tempName, "muPlus[%d]", obs); strcpy(reform->cname[k], tempName);
		reform->cname[k+1] = (cString) arr_alloc(NAMESIZE, char);
		sprintf(tempName, "muMinus[%d]", obs); strcpy(reform->cname[k+1], tempName);

		reform->matval[cnt] = -1;
		reform->matind[cnt++] = prob[0]->sp->mar + omega->numObs*prob[1]->sp->mar + obs;
		reform->matval[cnt] = 1;
		reform->matind[cnt++] = prob[0]->sp->mar + omega->numObs*prob[1]->sp->mar + obs;

		/* Add an additional constraint for moment matching */
		int idx = prob[0]->sp->mar + omega->numObs*prob[1]->sp->mar + obs;
		reform->rhsx[idx] = 0.0;
		reform->senx[idx] = 'E';
		reform->rname[idx] = (cString) arr_alloc(NAMESIZE, char);
		sprintf(tempName, "match[%d]", obs);
		strcpy(reform->rname[idx], tempName);
		k += 2;
	}

	/* Load the copy into CPLEX */
	reform->lp = setupProblem(reform->name, reform->type, reform->mac, reform->mar, reform->objsen, reform->objx, reform->rhsx, reform->senx, reform->matbeg,
			reform->matcnt,reform->matind, reform->matval, reform->bdl, reform->bdu, NULL, reform->cname, reform->rname, reform->ctype);
	if ( reform->lp == NULL ) {
		errMsg("Problem Setup", "setupMomentMatchReform", "failed to setup master problem in the solver",0);
		return NULL;
	}

#if defined(ALGO_CHECK)
	if ( writeProblem(reform->lp, "reformed.lp") ) {
		errMsg("solver", "setupMomentMatchReform", "failed to write master problem to file", 0);
		return NULL;
	}
#endif

	return reform;
}//END setupMomentMatchReform()

oneProblem *setupWassersteinInfReform (probType **prob, omegaType *omega) {
	oneProblem *reform = NULL;
	char tempName[NAMESIZE];

	reform = (oneProblem *) mem_malloc (sizeof(oneProblem));
	reform->type = PROB_LP;
	reform->objsen  = 1;

	reform->mar 	= prob[0]->sp->mar + omega->numObs*prob[1]->sp->mar;								/* number of rows */
	reform->mac 	= prob[0]->sp->mac + omega->numObs*prob[1]->sp->mac;          						/* number of columns + etas */
	reform->numInt 	= prob[0]->sp->numInt + omega->numObs*prob[1]->sp->numInt;    						/* number of integer variables in the problem  */
	reform->numnz 	= prob[0]->sp->numnz + omega->numObs*prob[1]->sp->numnz;      						/* number of non-zero elements in constraint matrix */
	reform->matsz 	= prob[0]->sp->matsz + omega->numObs*(prob[1]->sp->matsz + prob[1]->Cbar->cnt);	/* extended matrix size */
	reform->marsz 	= prob[0]->sp->marsz + omega->numObs*prob[1]->sp->marsz;      						/* extended row size */
	reform->macsz 	= prob[0]->sp->macsz + omega->numObs*prob[1]->sp->macsz;      						/* extended column size */

	reform->rstorsz = reform->cstorsz = 0;

	/* Allocate memory to the information whose type is cString */
	reform->name = (cString) arr_alloc(NAMESIZE, char);
	reform->senx = (cString) arr_alloc(reform->marsz,char);
	reform->ctype = (cString) arr_alloc(reform->macsz,char);
	reform->objname = (cString) arr_alloc(NAMESIZE,char);
	reform->cname = (cString*) arr_alloc(reform->macsz, cString);
	reform->rname = (cString *) arr_alloc(reform->marsz, cString);
	reform->cstore = NULL;
	reform->rstore = NULL;

	/* Allocate memory to the information whose type is dVector */
	reform->objx = (dVector) arr_alloc(reform->macsz, double);
	reform->rhsx = (dVector) arr_alloc(reform->marsz, double);
	reform->matval = (dVector) arr_alloc(reform->matsz, double);
	reform->bdl = (dVector) arr_alloc(reform->macsz, double);
	reform->bdu = (dVector) arr_alloc(reform->macsz, double);

	/* Allocate memory to the information whose type is iVector */
	reform->matbeg = (iVector) arr_alloc(reform->macsz, int);
	reform->matcnt = (iVector) arr_alloc(reform->macsz, int);
	reform->matind = (iVector) arr_alloc(reform->matsz, int);

	strcpy(reform->name, "WassersteinReform");           /* Copy problem name */
	strcpy(reform->objname, "WassersteinDistance");      /* Copy objective name */

	/* First-stage problem */
	int cnt = 0;
	/* Copy all the column information from the original first-stage problem */
	for (int j = 0; j < prob[0]->sp->mac; j++) {
		reform->objx[j] 	= prob[0]->sp->objx[j];					/* Copy objective function coefficients */
		reform->ctype[j] 	= prob[0]->sp->ctype[j];				/* Copy the decision variable type */
		reform->bdu[j] 		= prob[0]->sp->bdu[j];					/* Copy the upper bound and lower bound */
		reform->bdl[j] 		= prob[0]->sp->bdl[j];
		reform->matbeg[j] 	= cnt;									/* Copy the master sparse matrix beginning position of each column */
		reform->matcnt[j] = 0;
		for (int idx = prob[0]->sp->matbeg[j]; idx < prob[0]->sp->matbeg[j] + prob[0]->sp->matcnt[j]; idx++) {
			reform->matval[cnt] = prob[0]->sp->matval[idx];		/* Copy the non-zero coefficient */
			reform->matind[cnt] = prob[0]->sp->matind[idx];		/* Copy the row entry of the non-zero elements */
			reform->matcnt[j]++;
			cnt++;
		}

		/* Loop through the transfer matrix to introduce the first-stage variables for second-stage constraints for all observations */
		for (int idx = 0; idx < prob[1]->Cbar->cnt; idx++ ) {
			if ( prob[1]->Cbar->col[idx+1] == j + 1 ) {
				for ( int obs = 0; obs < omega->numObs; obs++ ) {
					reform->matval[cnt] = prob[1]->Cbar->val[idx+1];
					reform->matind[cnt] = prob[0]->sp->mar + obs*prob[1]->sp->mar + prob[1]->Cbar->row[idx+1]-1;
					reform->matcnt[j]++;
					cnt++;
				}
			}
		}

		reform->cname[j] = (cString) arr_alloc(NAMESIZE, char);
		strcpy(reform->cname[j], prob[0]->sp->cname[j]);
	}

	/* Copy all information concerning rows of first-stage */
	for (int r = 0; r < prob[0]->sp->mar; r++) {
		reform->rhsx[r]  = prob[0]->sp->rhsx[r]; 				/* Copy the right hand side value */
		reform->senx[r]  = prob[0]->sp->senx[r];				/* Copy the constraint sense */

		reform->rname[r] = (cString) arr_alloc(NAMESIZE, char);
		strcpy(reform->rname[r], prob[0]->sp->rname[r]);
	}

	/* Second-stage problem */
	/* Loop through all the observations to add the variables and constraints for each one. */
	for ( int obs = 0; obs < omega->numObs; obs++) {
		/* Copy all the column information from the original second-stage problem */
		int k = prob[0]->sp->mac + obs*prob[1]->sp->mac;
		for (int j = 0; j < prob[1]->sp->mac; j++) {
			reform->objx[k+j] 	= prob[1]->sp->objx[j]/(double) omega->numObs;					/* Copy objective function coefficients */
			reform->ctype[k+j] 	= prob[1]->sp->ctype[j];				/* Copy the decision variable type */
			reform->bdu[k+j] 	= prob[1]->sp->bdu[j];					/* Copy the upper bound and lower bound */
			reform->bdl[k+j] 	= prob[1]->sp->bdl[j];
			reform->matbeg[k+j] = cnt;									/* Copy the master sparse matrix beginning position of each column */
			reform->matcnt[k+j] = 0;
			for (int idx = prob[1]->sp->matbeg[j]; idx < prob[1]->sp->matbeg[j] + prob[1]->sp->matcnt[j]; idx++) {
				reform->matval[cnt] = prob[1]->sp->matval[idx];		/* Copy the non-zero coefficient */
				reform->matind[cnt] = prob[0]->sp->mar + obs*prob[1]->sp->mar + prob[1]->sp->matind[idx];		/* Copy the row entry of the non-zero elements */
				reform->matcnt[k+j]++;
				cnt++;
			}

			reform->cname[k+j] = (cString) arr_alloc(NAMESIZE, char);
			sprintf(tempName, "%s_%d", prob[1]->sp->cname[j], obs);
			strcpy(reform->cname[k+j], tempName);
		}

		/* Copy all information concerning rows of second-stage */
		k = prob[0]->sp->mar + obs*prob[1]->sp->mar;
		for (int r = 0; r < prob[1]->sp->mar; r++) {
			reform->rhsx[k+r]  = prob[1]->sp->rhsx[r]; 				/* Copy the right hand side value */
			reform->senx[k+r]  = prob[1]->sp->senx[r];				/* Copy the constraint sense */

			reform->rname[k+r] = (cString) arr_alloc(NAMESIZE, char);
			sprintf(tempName, "%s_%d", prob[1]->sp->rname[r], obs);
			strcpy(reform->rname[k+r], tempName);
		}
	}

	/* Update the right-hand side of constraints with random variables to rhs + theta */
	for ( int obs = 0; obs < omega->numObs; obs++ ) {
		for ( int rv = 0; rv < prob[1]->num->numRV; rv++ ) {
			int idx = prob[0]->sp->mar + obs*prob[1]->sp->mar + prob[1]->coord->rvRows[rv+1]-1;
			reform->rhsx[idx] = (omega->vals[obs][rv+1] + config.DRO_PARAM_2);
		}
	}

	/* Load the copy into CPLEX */
	reform->lp = setupProblem(reform->name, reform->type, reform->mac, reform->mar, reform->objsen, reform->objx, reform->rhsx, reform->senx, reform->matbeg,
			reform->matcnt,reform->matind, reform->matval, reform->bdl, reform->bdu, NULL, reform->cname, reform->rname, reform->ctype);
	if ( reform->lp == NULL ) {
		errMsg("Problem Setup", "setupWassersteinInfReform", "failed to setup master problem in the solver",0);
		return NULL;
	}

#if defined(ALGO_CHECK)
	if ( writeProblem(reform->lp, "reformed.lp") ) {
		errMsg("solver", "setupWassersteinInfReform", "failed to write master problem to file", 0);
		return NULL;
	}
#endif

	return reform;
}//END setupWassersteinInfReform()
