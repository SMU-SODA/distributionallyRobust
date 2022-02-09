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
int solveMaster(numType *num, sparseVector *dBar, cellType *cell, double lb) {
	int 	status;
	clock_t	tic;

	if ( config.ALGO_TYPE == SD ) {
		if( changeEtaCol(cell->master->lp, num->rows, num->cols, cell->omega->numObs, cell->cuts) ) {
			errMsg("algorithm", "solveQPMaster", "failed to change the eta column coefficients", 0);
			return 1;
		}
		cell->incumbChg = false;
	}



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
		cell->candidEst = vXvSparse(cell->candidX, dBar) + getPrimalPoint(cell->master->lp, num->cols);
	}
	else {
		cell->candidEst = getObjective(cell->master->lp, PROB_LP);
	}
	cell->gamma = cell->candidEst - cell->incumbEst;

	return 0;
}//END solveMaster()

/* This subroutine initializes the master problem by copying information from the decomposed prob[0](type: oneProblem) and adding a column for
 * theta for modified benders decomposition. */
oneProblem *newMaster(oneProblem *orig, double lb, int extraCols, double objCoeff) {
	oneProblem 	*master;
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
	master->cstorsz = orig->cstorsz + NAMESIZE*extraCols; 		/* memory size for storing column names */
	master->mac 	= orig->mac + extraCols;			/* number of columns + etas */
	master->macsz 	= orig->mac + extraCols;			/* extended column size */

	/* Allocate memory to the information whose type is cString */
	master->name = (cString) arr_alloc(NAMESIZE, char);
	master->senx = (cString) arr_alloc(master->marsz,char);
	master->ctype = (cString) arr_alloc(master->macsz,char);
	master->objname = (cString) arr_alloc(NAMESIZE,char);
	master->cname = (cString*) arr_alloc(master->macsz,cString);
	master->cstore = (cString) arr_alloc(master->cstorsz, char);
	if ( master->mar > 0 ) {
		master->rname = (cString *) arr_alloc(master->marsz,cString);
		master->rstore = (cString) arr_alloc(master->rstorsz, char);
	}
	else {
		master->rname = NULL; master->rstore = NULL;
	}

	/* Allocate memory to the information whose type is dVector */
	master->objx = (dVector) arr_alloc(master->macsz, double);
	master->rhsx = (dVector) arr_alloc(master->marsz, double);
	master->matval = (dVector) arr_alloc(master->matsz, double);
	master->bdl = (dVector) arr_alloc(master->macsz, double);
	master->bdu = (dVector) arr_alloc(master->macsz, double);

	/* Allocate memory to the information whose type is iVector */
	master->matbeg = (iVector) arr_alloc(master->macsz, int);
	master->matcnt = (iVector) arr_alloc(master->macsz, int);
	master->matind = (iVector) arr_alloc(master->matsz, int);

	strcpy(master->name, orig->name);           /* Copy problem name */
	strcpy(master->objname, orig->objname);     /* Copy objective name */

	/* Copy problem's column and row names, and calculate the pointers for master/copy row and column names. */
	int i = 0;
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
	int cnt = 0;
	for (int j = 0; j < orig->mac; j++) {
		master->objx[j] = orig->objx[j];		/* Copy objective function coefficients */
		master->ctype[j] = orig->ctype[j];		/* Copy the decision variable type */
		master->bdu[j] = orig->bdu[j];			/* Copy the upper bound and lower bound */
		master->bdl[j] = orig->bdl[j];
		master->cname[j] = orig->cname[j] + colOffset;	/* Copy column names, offset by length */
		master->matbeg[j] = cnt;				/* Copy the master sparse matrix beginning position of each column */
		master->matcnt[j] = orig->matcnt[j];	/* Copy the sparse matrix non-zero element count */
		master->ctype[j] = orig->ctype[j];		/* Loop through all non-zero elements in this column */
		for (int idx = orig->matbeg[j]; idx < orig->matbeg[j] + orig->matcnt[j]; idx++) {
			master->matval[cnt] = orig->matval[idx];	/* Copy the non-zero coefficient */
			master->matind[cnt] = orig->matind[idx];	/* Copy the row entry of the non-zero elements */
			cnt++;
		}
	}

	/* Copy all information concerning rows of master */
	for (int r = 0; r < orig->mar; r++) {
		master->rhsx[r] = orig->rhsx[r]; 	/* Copy the right hand side value */
		master->senx[r] = orig->senx[r];	/* Copy the constraint sense */
		master->rname[r] = orig->rname[r] + rowOffset; /* Copy row names, offset by length */
	}

	/* Initialize information for the extra columns in the new master. */
	char tempName[NAMESIZE];
	colOffset = orig->cstorsz;
	for ( int n = 0; n < extraCols; n++ ) {
		sprintf(tempName, "eta_%d", n);
		strcpy(master->cstore + colOffset, tempName);
		master->cname[orig->mac+n] = master->cstore + colOffset;
		master->objx[orig->mac+n] = (n == 0) ? objCoeff:1.0/(double) (extraCols-1);
		master->ctype[orig->mac+n] = 'C';
		master->bdu[orig->mac+n] = INFBOUND;
		master->bdl[orig->mac+n] = lb;
		master->matbeg[orig->mac+n] = orig->numnz;
		master->matcnt[orig->mac+n] = 0;

		colOffset += strlen(tempName) + 1;
	}

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

int addCut2Master(cellType *cell, cutsType *cuts, oneCut *cut, int lenX, int obsID) {
	iVector 	indices;
	int 	cnt;
	static int cummCutNum = 0;

	/* If it is optimality cut being added, check to see if there is room for the candidate cut, else drop a cut */
	if (cuts->cnt == cell->maxCuts) {
		/* make room for the latest cut */
		if( reduceCuts() < 0 ) {
			errMsg("algorithm", "addCut2Master", "failed to add reduce cuts to make room for candidate cut", 0);
			return -1;
		}
	}

	cut->alphaIncumb = cut->alpha;
	if ( config.MASTER_TYPE == PROB_QP )
		cut->alphaIncumb -= vXv(cut->beta, cell->incumbX, NULL, lenX);

	if (!(indices = arr_alloc(lenX + 1, int)))
		errMsg("Allocation", "addcut2Master", "fail to allocate memory to coefficients of beta",0);
	for (cnt = 1; cnt <= lenX; cnt++)
		indices[cnt] = cnt - 1;
	indices[0] = lenX;

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

/* This function performs the updates on all the coefficients of eta in the master problem constraint matrix.  During every iteration,
 * each of the coefficients on eta are increased, so that the effect of the cut on the objective function is decreased. */
int changeEtaCol(LPptr lp, int numRows, int numCols, int currObs, cutsType *cuts) {
	double	coef[1];
	int 	c;

	for (c = 0; c < cuts->cnt; c++){
		/* Currently both incumbent and candidate cuts are treated similarly, and sunk as iterations proceed */
		coef[0] = (double) (currObs) / (double) cuts->vals[c]->numObs;         // coefficient k/j of eta column

		if ( changeCol(lp, numCols, coef, cuts->vals[c]->rowNum, cuts->vals[c]->rowNum+1) ) {
			errMsg("solver", "changeEtaCol", "failed to change eta column in the stage problem", 0);
			return 1;
		}
	}

	return 0;
}//END changeEtaCol()

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

int updateRHS(LPptr lp, cutsType *cuts, int numIter, double lb) {
	int 	cnt;
	dVector	rhs;
	iVector	indices;

	if (!(rhs = arr_alloc(cuts->cnt, double)))
		errMsg("allocation", "updateRHS", "rhs", 0);
	if (!(indices = arr_alloc(cuts->cnt, int)))
		errMsg("allocation", "updateRHS", "indices", 0);

	for (cnt = 0; cnt < cuts->cnt; cnt++) {
		rhs[cnt] = cuts->vals[cnt]->alphaIncumb + ((double) numIter / (double) cuts->vals[cnt]->numObs - 1) * lb;
		indices[cnt] = cuts->vals[cnt]->rowNum;
	}

	/* Now we change the right-hand of the master problem. */
	if ( changeRHS(lp, cuts->cnt, indices, rhs) ) {
		errMsg("solver", "updateRHS", "failed to change the right-hand side in the solver", 0);
		return 1;
	}

	mem_free(rhs);
	mem_free(indices);

	return 0;
}//END updateRHS

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
