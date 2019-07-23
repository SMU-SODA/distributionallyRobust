/*
 * distributeSep.c
 *
 *  Created on: Apr 9, 2019
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "samplingDRO.h"

oneProblem *newDistributeSep(int numMoments) {
	oneProblem *sep;

	if ( (sep = (oneProblem *) malloc(sizeof(oneProblem))) )
		errMsg("memory allocation", "newDistributeSep", "failed to allocate memory to oneProblem", 0);

	/* Assign values to elements of oneProblem */
	sep->type 	 = PROB_LP;   		/* type of problem: LP, QP, MIP or MIQP */
	sep->objsen  = SOLVER_MIN;      /* sense of the objective: 1 for minimization and -1 for maximization */
	sep->mar 	 = 0;               /* number of rows */
	sep->numInt  = 0;			    /* number of integer variables in the problem  */
	sep->numnz 	 = 0;               /* number of non-zero elements in constraint matrix */
	sep->matsz 	 = 0;               /* extended matrix size */
	sep->marsz 	 = 0;               /* extended row size */
	sep->rstorsz = 0;               /* memory size for storing row names */
	sep->mac 	 = 0;  				/* number of columns + etas */
	sep->macsz 	 = 0;				/* extended column size */
	sep->cstorsz = 0;  				/* memory size for storing column names */

	/* Allocate memory to the information whose type is cString */
	if (!(sep->name = (cString) arr_alloc(NAMESIZE, char)))
		errMsg("Allocation", "newDistributeSep", "Fail to allocate memory to master->name",0);


	return sep;
}//END newDistributeSep()
