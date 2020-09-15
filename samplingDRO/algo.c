/*
 * algo.c
 *
 *  Created on: Mar 22, 2019
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */
#include "samplingDRO.h"

extern configType config;
extern cString outputDir;

int algo(oneProblem *orig, timeType *tim, stocType *stoc, cString inputDir, cString probName) {
	probType **prob = NULL;
	cellType *cell = NULL;
	dVector  meanSol;
	clock_t	 tic;
	int		 rep;
	FILE 	*detailedResults = NULL, *summary = NULL, *incumbFile = NULL;

	/* Open solver */
	openSolver();

	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell, &meanSol) )
		goto TERMINATE;

	detailedResults = openFile(outputDir, "detailedResults.csv", "w");
	incumbFile = openFile(outputDir, "incumb.dat", "w");
	summary = openFile(outputDir, "summary.dat", "w");

	printDecomposeSummary(summary, probName, tim, prob);
	printDecomposeSummary(stdout, probName, tim, prob);

	for ( rep = 0; rep < config.MULTIPLE_REP; rep++ ) {
		fprintf(stdout, "\n====================================================================================================================================\n");
		fprintf(stdout, "Replication-%d\n", rep+1);

		/* setup the seed to be used in the current iteration */
		config.RUN_SEED[0]  = config.RUN_SEED[rep+1];
		config.EVAL_SEED[0] = config.EVAL_SEED[rep+1];

		if ( rep > 0 ) {
			/* clean up the cell for the next replication */
			if ( cleanCellType(cell, prob[0], meanSol) ) {
				errMsg("algorithm", "benders", "failed to solve the cells using MASP algorithm", 0);
				goto TERMINATE;
			}
		}

		/* Update omega structure */
		if ( config.SAMPLING_TYPE == 1 ) {
			cell->omega->cnt = config.MAX_OBS;
			setupSAA(stoc, NULL, &config.RUN_SEED[0], &cell->omega->vals, &cell->omega->probs, &cell->omega->cnt,
					config.TOLERANCE);
			for ( int m = 0; m < cell->omega->cnt; m++ )
				for ( int n = 1; n <= stoc->numOmega; n++ )
					cell->omega->vals[m][n] -= stoc->mean[n-1];
		}

		/* Setup the distribution separation problem */
		cell->sep = newDistSepProb(stoc, cell->omega);

		tic = clock();
		/* Use two-stage algorithm to solve the problem */
		if ( solveFixedDROCell(stoc, prob, cell) ) {
			errMsg("algorithm", "benders", "failed to solve the cells using MASP algorithm", 0);
			goto TERMINATE;
		}
		cell->time->repTime = ((double) (clock() - tic))/CLOCKS_PER_SEC;

		/* Write solution statistics for optimization process */
		printOptimizationSummary(cell);
		writeOptimizationStatistics(detailedResults, incumbFile, prob, cell, rep);

		/* evaluating the optimal solution*/
		if (config.EVAL_FLAG == 1) {
			dVector evalX = (config.MASTER_TYPE == PROB_QP || config.SAMPLING_TYPE == 2) ? cell->incumbX : cell->candidX;
			evaluate(detailedResults, stoc, prob, cell, evalX);
		}
		else
			fprintf(detailedResults,"\n");
	}

	printf("\n\nSuccessfully completed the DR-SD algorithm.\n");

	/* release structures and close solver environment */
	if (cell) freeCellType(cell);
	if (prob) freeProbType(prob,2);
	if (orig) freeOneProblem(orig);
	if (tim) freeTimeType(tim);
	if (stoc)freeStocType(stoc);
	if (meanSol) mem_free(meanSol);
	closeSolver();
	return 0;

	TERMINATE:
	if (cell) freeCellType(cell);
	if (prob) freeProbType(prob,2);
	if (orig) freeOneProblem(orig);
	if (tim) freeTimeType(tim);
	if (stoc)freeStocType(stoc);
	if (meanSol) mem_free(meanSol);
	closeSolver();
	return 1;
}//END algo()

int solveFixedDROCell(stocType *stoc, probType **prob, cellType *cell) {
	int candidCut;
	clock_t tic;

	/* Main loop of the algorithm */
	while ( cell->k < config.MAX_ITER ) {
		tic = clock();
		cell->k++;

		/******* 1. Solve the subproblem with candidate solution, form and update the candidate cut *******/
		if ( (candidCut = formOptCut(prob[1], cell, cell->candidX)) < 0 ) {
			errMsg("algorithm", "solveCell", "failed to add candidate cut", 0);
			return 1;
		}

		/******* 2. Check for optimality of the current solution *******/
		if (optimal(prob[0], cell)) {
			cell->incumbEst = cell->candidEst;
			break;
		}

		/******* 3. Solve the master problem to obtain the new candidate solution *******/
		if ( solveMaster(prob[0]->num, prob[0]->dBar, cell) ) {
			errMsg("algorithm", "solveCell", "failed to solve master problem", 0);
			return 1;
		}

		/* Update the solution times */
		cell->time->masterAccumTime += cell->time->masterIter; cell->time->subprobAccumTime += cell->time->subprobIter;
		cell->time->argmaxAccumTime += cell->time->argmaxIter; cell->time->optTestAccumTime += cell->time->optTestIter;
		cell->time->masterIter = cell->time->subprobIter = cell->time->optTestIter = 0.0;
		cell->time->argmaxIter = cell->time->optTestIter = 0.0;
		cell->time->iterTime = ((double) clock() - tic)/CLOCKS_PER_SEC; cell->time->iterAccumTime += cell->time->iterTime;
	}//END main loop

	return 0;
}//END solveFixedDROCell()

int solveSeqDROCell(stocType *stoc, probType **prob, cellType *cell) {
	clock_t tic;
	int candidCut;

	/* Main loop of the algorithm */
	while ( cell->k < config.MAX_ITER ) {
		tic = clock();
		cell->k++;

		/******* 1. Check for optimality of the current solution *******/
		if ( config.MASTER_TYPE == PROB_QP ) {
			if ( optimal(prob[0], cell) ) {
				break;
			}
		}

		/******* 2. Solve the subproblem with candidate solution, form and update the candidate cut *******/
		if ( (candidCut = formOptCut(prob[1], cell, cell->candidX)) < 0 ) {
			errMsg("algorithm", "solveCell", "failed to add candidate cut", 0);
			return 1;
		}


	}//END main loop

	return 0;
}//END solveSeqDROPCell()

void writeOptimizationStatistics(FILE *soln, FILE *incumb, probType **prob, cellType *cell, int rep) {

	/* Print header for the first replication*/
	if ( rep == 0)
		fprintf(soln, "Replication\tIterations\tOpt estimate\tTotal time\tMaster time\t Subproblem time\t Optimality time\tArgmax time\t Reduce time\t"
				"Eval Estimate\tError\tCI-L\tCI-U\tOutcomes\n");

	fprintf(soln, "%d\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf", rep+1, cell->k, cell->incumbEst,cell->time->repTime, cell->time->masterAccumTime,
			cell->time->subprobAccumTime, cell->time->optTestAccumTime, cell->time->argmaxAccumTime, cell->time->reduceTime);

	if ( incumb != NULL ) {
		if ( config.MASTER_TYPE == PROB_QP )
			printVector(cell->incumbX, prob[0]->num->cols, incumb);
		else
			printVector(cell->candidX, prob[0]->num->cols, incumb);
	}

}//END writeOptimizationStatistics()

/* Prints summary statistics to the console output. */
void printOptimizationSummary(cellType *cell) {

	fprintf(stdout, "\n------------------------------------------------------------ Optimization ---------------------------------------------------------\n");
	if ( config.MASTER_TYPE == PROB_QP )
		fprintf(stdout, "Algorithm                          : Regularized Decomposition-based DRO \n");
	else
		fprintf(stdout, "Algorithm                          : Decomposition-based DRO \n");
	fprintf(stdout, "Number of iterations               : %d", cell->k);
	if ( cell->k == config.MAX_ITER)
		fprintf(stdout, "*\n");
	else
		fprintf(stdout, "\n");
	fprintf(stdout, "Optimization value                 : %f\n", cell->incumbEst);
	fprintf(stdout, "Total time                         : %f\n", cell->time->repTime);
	fprintf(stdout, "Total time to solve master         : %f\n", cell->time->masterAccumTime);
	fprintf(stdout, "Total time to solve subproblems    : %f\n", cell->time->subprobAccumTime);
	fprintf(stdout, "Total time to verify optimality    : %f\n", cell->time->optTestAccumTime);

}//END writeOptimizationSummary()
