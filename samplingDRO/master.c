/*
 * master.c
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

/* In this his function the master problem is solved after the newest cut is added to master problem, the incumbent cut is updated if necessary.
 * Here the coefficients on all the cuts are updated, and finally master problem is solved. */
int solveMaster(numType *num, sparseVector *dBar, cellType *cell) {
	int 	status;
	clock_t	tic;

#if defined(ALGO_CHECK)
	writeProblem(cell->master->lp,"cellMaster.lp");
#endif

	tic = clock();
	/* solve the master problem */
	if ( solveProblem(cell->master->lp, cell->master->name, config.MASTER_TYPE, &status) ) {
		writeProblem(cell->master->lp, "error.lp");
		errMsg("algorithm", "solveMaster", "failed to solve the master problem", 0);
		return 1;
	}
	cell->time->masterIter = ((double) (clock() - tic))/CLOCKS_PER_SEC;

	/* increment the number of problems solved during algorithm */
	cell->LPcnt++;

	/* Get the most recent optimal solution to master program */
	if ( getPrimal(cell->master->lp, cell->candidX, num->cols) ) {
		errMsg("algorithm", "solveMaster", "failed to obtain the primal solution for master", 0);
		return 1;
	}

	if ( cell->master->type == PROB_QP ) {
		/* Get the dual solution too */
		if ( getDual(cell->master->lp, cell->piM, cell->master->mar) ) {
			errMsg("solver", "solveQPMaster", "failed to obtain dual solutions to master", 0);
			return 1;
		}

		addVectors(cell->candidX, cell->incumbX, NULL, num->cols);
		cell->candidEst = vXvSparse(cell->candidX, dBar);
		cell->candidEst += getPrimalPoint(cell->master->lp, num->cols);
		cell->gamma = cell->candidEst - cell->incumbEst;
	}
	else
		cell->candidEst = getObjective(cell->master->lp, PROB_LP);

	return 0;
}//END solveMaster()

/* This subroutine initializes the master problem by copying information from the decomposed prob[0](type: oneProblem) and adding a column for
 * theta for modified benders decomposition. */
oneProblem *newMaster(oneProblem *orig, double lb, omegaType *omega) {
	oneProblem 	*master;
	int         r, i, j, idx, cnt;
	long        colOffset, rowOffset;
	char        *q;

	if (!(master = (oneProblem *) mem_malloc (sizeof(oneProblem))))
		errMsg("Memory allocation", "new_master", "Failed to allocate memory to mcell->sp", 0);

	/* -+-+-+-+-+-+-+-+-+-+-+-+-+-+- Allocating memory to master -+-+-+-+-+-+-+-+-+-+-+-+-+-+- */
	master->type 	= config.MASTER_TYPE;               /* type of problem: LP, QP, MIP or MIQP */
	master->objsen 	= orig->objsen;                 	/* sense of the objective: 1 for minimization and -1 for maximization */
	master->mar 	= orig->mar;                       	/* number of rows */
	master->numInt 	= orig->numInt;                 	/* number of integer variables in the problem  */
	master->numnz 	= orig->numnz;                   	/* number of non-zero elements in constraint matrix */
	master->matsz 	= orig->matsz;                   	/* extended matrix size */
	master->marsz 	= orig->marsz;                   	/* extended row size */
	master->rstorsz = orig->rstorsz;               		/* memory size for storing row names */
	master->mac 	= orig->mac+1;           			/* number of columns + etas */
	master->macsz 	= orig->macsz + 1;       			/* extended column size */
	master->cstorsz = orig->cstorsz + NAMESIZE;    		/* memory size for storing column names */

	/* Allocate memory to the information whose type is cString */
	if (!(master->name = (cString) arr_alloc(NAMESIZE, char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->name",0);
	if (!(master->senx = (cString) arr_alloc(master->marsz,char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->senx",0);
	if (!(master->ctype = (cString) arr_alloc(master->macsz,char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->ctype",0);
	if (!(master->objname = (cString) arr_alloc(NAMESIZE,char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->objname",0);
	if (!(master->cname = (cString*) arr_alloc(master->macsz,cString)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->cname",0);
	if (!(master->cstore = (cString) arr_alloc(master->cstorsz, char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->cstore",0);
	if ( master->mar > 0 ) {
		if (!(master->rname = (cString *) arr_alloc(master->marsz,cString)))
			errMsg("Allocation", "new_master", "Fail to allocate memory to master->rname",0);
		if (!(master->rstore = (cString) arr_alloc(master->rstorsz, char)))
			errMsg("Allocation", "new_master", "Fail to allocate memory to master->rstore",0);
	}
	else {
		master->rname = NULL; master->rstore = NULL;
	}

	/* Allocate memory to the information whose type is dVector */
	if (!(master->objx = (dVector) arr_alloc(master->macsz, double)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->objx",0);
	if (!(master->rhsx = (dVector) arr_alloc(master->marsz, double)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->rhsx",0);
	if (!(master->matval = (dVector) arr_alloc(master->matsz, double)))
		errMsg("allocation", "new_master", "master->matval",0);
	if (!(master->bdl = (dVector) arr_alloc(master->macsz, double)))
		errMsg("allocation", "new_master", "master->bdl",0);
	if (!(master->bdu = (dVector) arr_alloc(master->macsz, double)))
		errMsg("allocation", "new_master", "master->bdu",0);

	/* Allocate memory to the information whose type is iVector */
	if (!(master->matbeg = (iVector) arr_alloc(master->macsz, int)))
		errMsg("allocation", "new_master", "master->matbeg",0);
	if (!(master->matcnt = (iVector) arr_alloc(master->macsz, int)))
		errMsg("allocation", "new_master", "master->matcnt",0);
	if (!(master->matind = (iVector) arr_alloc(master->matsz, int)))
		errMsg("allocation", "new_master", "master->matind",0);

	strcpy(master->name, orig->name);           /* Copy problem name */
	strcpy(master->objname, orig->objname);     /* Copy objective name */

	/* Copy problem's column and row names, and calculate the pointers for master/copy row and column names. */
	i = 0;
	for (q = orig->cname[0]; q < orig->cname[0] + orig->cstorsz; q++)
		master->cstore[i++] = *q;
	colOffset = master->cstore - orig->cname[0];

	if ( master->mar > 0 ) {
		i = 0;
		for (q = orig->rname[0]; q < orig->rname[0] + orig->rstorsz; q++)
			master->rstore[i++] = *q;
		rowOffset = master->rstore - orig->rname[0];
	}

	/* Copy all the column information from the original master problem */
	cnt = 0;
	for (j = 0; j < orig->mac; j++) {
		/* Copy objective function coefficients */
		master->objx[j] = orig->objx[j];
		/* Copy the decision variable type */
		master->ctype[j] = orig->ctype[j];
		/* Copy the upper bound and lower bound */
		master->bdu[j] = orig->bdu[j];
		master->bdl[j] = orig->bdl[j];
		/* Copy column names, offset by length */
		master->cname[j] = orig->cname[j] + colOffset;
		/* Copy the master sparse matrix beginning position of each column */
		master->matbeg[j] = cnt;
		/* Copy the sparse matrix non-zero element count */
		master->matcnt[j] = orig->matcnt[j];
		master->ctype[j] = orig->ctype[j];
		/* Loop through all non-zero elements in this column */
		for (idx = orig->matbeg[j]; idx < orig->matbeg[j] + orig->matcnt[j]; idx++) {
			/* Copy the non-zero coefficient */
			master->matval[cnt] = orig->matval[idx];
			/* Copy the row entry of the non-zero elements */
			master->matind[cnt] = orig->matind[idx];
			cnt++;
		}
	}

	/* Copy all information concerning rows of master */
	for (r = 0; r < orig->mar; r++) {
		/* Copy the right hand side value */
		master->rhsx[r] = orig->rhsx[r];
		/* Copy the constraint sense */
		master->senx[r] = orig->senx[r];
		/* Copy row names, offset by length */
		master->rname[r] = orig->rname[r] + rowOffset;
	}

	/* Initialize information for the extra column in the new master. */
	colOffset = orig->cstorsz;
	strcpy(master->cstore + orig->cstorsz, "eta");
	master->cname[orig->mac] = master->cstore + colOffset;
	master->objx[orig->mac] = 1.0;			// orig->mac is the last column in the original master
	master->ctype[orig->mac] = 'C';
	master->bdu[orig->mac] = INFBOUND;
	master->bdl[orig->mac] = lb;
	master->matbeg[orig->mac] = orig->numnz;	// Beginning point in matval/matind in eta columns. every eta column begins at the same address
	master->matcnt[orig->mac] = 0;               // Only optimality cuts has eta

	/* Load the copy into CPLEX */
	master->lp = setupProblem(master->name, master->type, master->mac, master->mar, master->objsen, master->objx, master->rhsx, master->senx, master->matbeg, master->matcnt,master->matind, master->matval, master->bdl, master->bdu, NULL, master->cname, master->rname, master->ctype);
	if ( master->lp == NULL ) {
		errMsg("Problem Setup", "new_master", "failed to setup master problem in the solver",0);
		return NULL;
	}

#if defined(ALGO_CHECK)
	if ( writeProblem(master->lp, "newMaster.lp") ) {
		errMsg("solver", "newMaster", "failed to write master problem to file", 0);
		return NULL;
	}
#endif

	return master;

}//END newMaster()

int checkImprovement(probType *prob, cellType *cell, int candidCut) {

	/* Calculate height at both the candidate and incumbent solution need to be computed. */
	cell->candidEst = vXvSparse(cell->candidX, prob->dBar) +
			maxCutHeight(cell->cuts, cell->candidX, prob->num->cols, 0);

#if defined(ALGO_CHECK)
	printf("Incumbent estimate = %lf; Candidate estimate = %lf\n", cell->incumbEst, cell->candidEst);
#endif

	/* If we see considerable improvement, then change the incumbent */
	if ( (cell->candidEst - cell->incumbEst) <= config.R1*cell->gamma ) {
		/* when we find an improvement, then we need to replace the incumbent x with candidate x */
		if ( replaceIncumbent(prob, cell) ) {
			errMsg("algorithm", "checkImprovement", "failed to replace incumbent solution with candidate", 0);
			return 1;
		}

		cell->cuts->vals[cell->iCutIdx[0]]->isIncumb = false;
		cell->iCutIdx[0] = candidCut;
		cell->cuts->vals[candidCut]->isIncumb = true;

		cell->incumbChg = false;
		printf("+"); fflush(stdout);
	}

	return 0;
}//END checkImprovement()

int replaceIncumbent(probType *prob, cellType *cell) {

	/* replace the incumbent solution with the candidate solution */
	copyVector(cell->candidX, cell->incumbX, prob->num->cols, 1);
	cell->incumbEst = cell->candidEst;

	/* update the right-hand side and the bounds with new incumbent solution */
	if ( constructQP(prob, cell, cell->incumbX, cell->quadScalar) ) {
		errMsg("algorithm", "replaceIncumbent", "failed to change the right-hand side after incumbent change", 0);
		return 1;
	}

	/* update the candidate cut as the new incumbent cut */
	cell->incumbChg = true;

	/* Since incumbent solution is now replaced by a candidate, we assume it is feasible now */
	cell->infeasIncumb = false;

	return 0;
}//END replaceIncumbent()

int constructQP(probType *prob, cellType *cell, dVector incumbX, double quadScalar) {
	int status;

	status = changeQPproximal(cell->master->lp, prob->num->cols, quadScalar, 0);
	if ( status ) {
		errMsg("algorithm", "algoIntSD", "failed to change the proximal term", 0);
		return 1;
	}
	status = changeQPrhs(prob, cell, incumbX);
	if ( status ) {
		errMsg("algorithm", "algoIntSD", "failed to change the right-hand side to convert the problem into QP", 0);
		return 1;
	}
	status = changeQPbds(cell->master->lp, prob->num->cols, prob->sp->bdl, prob->sp->bdu, incumbX);
	if ( status ) {
		errMsg("algorithm", "algoIntSD", "failed to change the bounds to convert the problem into QP", 0);
		return 1;
	}

	return status;
}//END constructQP()

/* Construct the Q diagonal matrix and copy it for quadratic problem. */
int changeQPproximal(LPptr lp, int numCols, double sigma, int numEta) {
	int    n;
	dVector qsepvec;

	if (!(qsepvec = arr_alloc(numCols+numEta+1, double)))
		errMsg("Allocation", "constructQP", "qsepvec",0);

	/* Construct Q matrix, which is simply a diagonal matrix. */
	for (n = 0; n < numCols; n++)
		qsepvec[n] = 0.5 * sigma;
	for ( n = numCols; n < numCols+numEta+1; n++ )
		qsepvec[n] = 0.0;

	/* Now copy the Q matrix for QP problem. */
	if ( copyQPseparable(lp, qsepvec) ) {
		errMsg("solver", "constructQP", "failed to copy Q matrix", 0);
		return 1;
	}

	mem_free(qsepvec);
	return 0;
}//END changeQPproximal()

/* In the regularized QP method, we need to change the rhs of x to d. The
 * 		 A * x 			= b
 * 		 eta + beta * x >= alpha
 * Since x = xbar + d, the corresponding changes will therefore be:
 * 		 A * d = b - A * xbar
 * 		 eta + beta * d >= alpha - beta * xbar
 * But as long as the incumbent sulotion does not change, b - A * xbar and alpha - beta * xbar (for the existing cuts) won't change. So we only need
 * to change it when the incumbent changes.
 *
 * On the other hand, in each iteration, a new cut will be added (and/or some cuts may be dropped) and therefore we need to shift the rhs of the
 * added cut from _alpha_ to _alpha - beta * xbar_, which has taken care of in the routine addCut() in cuts.c. We do not need to worry about the shift
 * of rhs for the dropped cuts.
 * This function performs the change of rhs when the incumbent changes, as described above. */
int changeQPrhs(probType *prob, cellType *cell, dVector xk) {
	int 	status = 0, cnt;
	dVector 	rhs;
	iVector 	indices;

	if (!(rhs =(dVector) arr_alloc(prob->num->rows+cell->cuts->cnt+1, double)))
		errMsg("Allocation", "changeRhs", "rhs",0);
	if (!(indices =(iVector) arr_alloc(prob->num->rows+cell->cuts->cnt, int)))
		errMsg("Allocation", "changeRhs", "indices",0);
	/* Be careful with the one_norm!! In the CxX() routine, it assumes the 0th element is reserved for the 1_norm, in the returned dVector, the T sparse
	 dVector, and the x dVector. */
	for (cnt = 0; cnt < prob->num->rows; cnt++) {
		rhs[cnt + 1] = prob->sp->rhsx[cnt];
		indices[cnt] = cnt;
	}

	/* b - A * xbar */
	rhs = MSparsexvSub(prob->Dbar, xk, rhs);

	/*** new rhs = alpha - beta * xbar (benders cuts)***/
	for (cnt = 0; cnt < cell->cuts->cnt; cnt++) {
		rhs[prob->num->rows+cnt+1] = cell->cuts->vals[cnt]->alpha - vXv(cell->cuts->vals[cnt]->beta, xk, NULL, prob->sp->mac);
		indices[prob->num->rows+cnt] = cell->cuts->vals[cnt]->rowNum;

		cell->cuts->vals[cnt]->alphaIncumb = rhs[prob->num->rows+cnt+1];
	}

	/* Now we change the right-hand of the master problem. */
	status = changeRHS(cell->master->lp, prob->num->rows + cell->cuts->cnt, indices, rhs+1);
	if (status)	{
		errMsg("solver", "changeQPrhs", "failed to change the right-hand side in the solver", 0);
		return 1;
	}

	mem_free(rhs);
	mem_free(indices);
	return 0;
}//END changeQPrhs()

/* This function changes the (lower) bounds of the variables, while changing from x to d. The lower bounds of d varibles are -xbar
 * (incumbent solution). */
int changeQPbds(LPptr lp, int numCols, dVector bdl, dVector bdu, dVector xk) {
	int 	status = 0, cnt;
	dVector	lbounds, ubounds;
	iVector	lindices, uindices;
	char 	*llu, *ulu;

	if (!(lbounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "changeBounds", "lbounds",0);
	if (!(lindices = arr_alloc(numCols, int)))
		errMsg("Allocation", "change_bounds", "lindices",0);
	if (!(llu = arr_alloc(numCols, char)))
		errMsg("Allocation", "changeBounds", "llu",0);

	if (!(ubounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "change_bounds", "ubounds",0);
	if (!(uindices = arr_alloc(numCols, int)))
		errMsg("Allocation", "changeBounds", "uindices",0);
	if (!(ulu = arr_alloc(numCols, char)))
		errMsg("Allocation", "changeBounds", "ulu",0);

	/* Change the Upper Bound */
	for (cnt = 0; cnt < numCols; cnt++) {
		ubounds[cnt] = bdu[cnt] - xk[cnt + 1];
		uindices[cnt] = cnt;
		ulu[cnt] = 'U';
	}

	status = changeBDS(lp, numCols, uindices, ulu, ubounds);
	if (status) {
		errMsg("algorithm", "changeQP", "failed to change the upper bound in the solver", 0);
		return 1;
	}

	/* Change the Lower Bound */
	for (cnt = 0; cnt < numCols; cnt++) {
		lbounds[cnt] = bdl[cnt] - xk[cnt + 1];
		lindices[cnt] = cnt;
		llu[cnt] = 'L';
	}

	status = changeBDS(lp, numCols, lindices, llu, lbounds);
	if (status) {
		errMsg("algorithm", "changeQP", "failed to change the lower bound in the solver", 0);
		return 1;
	}

	mem_free(lbounds); mem_free(lindices); mem_free(llu);
	mem_free(ubounds); mem_free(uindices); mem_free(ulu);

	return 0;
}//END changeQPbds()

int addCut2Master(cellType *cell, cutsType *cuts, oneCut *cut, int lenX, int obsID) {
	iVector 	indices;
	int 	cnt;
	static int cummCutNum = 0;

	/* If it is optimality cut being added, check to see if there is room for the candidate cut, else drop a cut */
	if (cuts->cnt == cell->maxCuts) {
		/* make room for the latest cut */
		if( reduceCut(cell->master, cuts, cell->candidX, cell->piM, lenX, cell->iCutIdx,
				cell->omega, obsID) < 0 ) {
			errMsg("algorithm", "addCut2Master", "failed to add reduce cuts to make room for candidate cut", 0);
			return -1;
		}
	}

	if ( config.MASTER_TYPE == PROB_QP )
		cut->alphaIncumb = cut->alpha - vXv(cut->beta, cell->incumbX, NULL, lenX);

	if (!(indices = arr_alloc(lenX + 1, int)))
		errMsg("Allocation", "addcut2Master", "fail to allocate memory to coefficients of beta",0);
	for (cnt = 1; cnt <= lenX; cnt++)
		indices[cnt] = cnt - 1;
	indices[0] = lenX+cut->omegaID;

	/* Add the cut to the cell cuts structure and assign a row number. */
	cuts->vals[cuts->cnt] = cut;
	cut->rowNum = cell->master->mar++;

	/* Set up the cut name */
	sprintf(cut->name, "cut_%04d", cummCutNum++);

	/* Add the row in the solver */
	if ( addRow(cell->master->lp, lenX + 1, cut->alphaIncumb, GE, 0, indices, cut->beta, cut->name) ) {
		errMsg("solver", "addcut2Master", "failed to add new row to problem in solver", 0);
		return -1;
	}

	mem_free(indices);
	return cuts->cnt++;
}//END addCuts2Master()
