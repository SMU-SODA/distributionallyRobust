/*
 * setup.c
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

int readConfig() {
	FILE 	*fptr;
	char	line[2*BLOCKSIZE], comment[2*BLOCKSIZE];
	int 	status, maxReps = 30, numReps = 0, numEvals = 0;

	fptr = fopen("config.dro", "r");
	if ( fptr == NULL ) {
		errMsg("read", "readConfig", "failed to open configuration file", 0);
		return 1;
	}

	config.RUN_SEED = (long long *) arr_alloc(maxReps+1, long long);
	config.EVAL_SEED = (long long *) arr_alloc(maxReps+1, long long);
	config.MAX_OBS = 0;

	while ((status = (fscanf(fptr, "%s", line) != EOF))) {
		if (!(strcmp(line, "RUN_SEED"))) {
			fscanf(fptr, "%lld", &config.RUN_SEED[numReps+1]);
			numReps++;
			if ( numReps > maxReps + 1 ) {
				maxReps *= 2;
				config.RUN_SEED = (long long *) mem_realloc(config.RUN_SEED, (maxReps+1)*sizeof(long long));
			}
		}
		else if (!(strcmp(line, "TOLERANCE")))
			fscanf(fptr, "%lf", &config.TOLERANCE);
		else if (!(strcmp(line, "MIN_ITER")))
			fscanf(fptr, "%d", &config.MIN_ITER);
		else if (!(strcmp(line, "MAX_ITER")))
			fscanf(fptr, "%d", &config.MAX_ITER);
		else if (!(strcmp(line, "MASTER_TYPE")))
			fscanf(fptr, "%d", &config.MASTER_TYPE);
		else if (!(strcmp(line, "CUT_MULT")))
			fscanf(fptr, "%d", &config.CUT_MULT);
		else if (!(strcmp(line, "R1")))
			fscanf(fptr, "%lf", &config.R1);
		else if (!(strcmp(line, "R2")))
			fscanf(fptr, "%lf", &config.R2);
		else if (!(strcmp(line, "R3")))
			fscanf(fptr, "%lf", &config.R3);
		else if (!(strcmp(line, "EPSILON")))
			fscanf(fptr, "%lf", &config.EPSILON);

		else if (!(strcmp(line, "CUT_MULT")))
			fscanf(fptr, "%d", &config.CUT_MULT);
		else if (!(strcmp(line, "TAU")))
			fscanf(fptr, "%d", &config.TAU);

		else if (!(strcmp(line, "MIN_QUAD_SCALAR")))
			fscanf(fptr, "%lf", &config.MIN_QUAD_SCALAR);
		else if (!(strcmp(line, "MAX_QUAD_SCALAR")))
			fscanf(fptr, "%lf", &config.MAX_QUAD_SCALAR);

		else if (!(strcmp(line, "MULTIPLE_REP")))
			fscanf(fptr, "%d", &config.MULTIPLE_REP);

		else if (!(strcmp(line, "SAA")))
			fscanf(fptr, "%d", &config.SAA);
		else if (!(strcmp(line, "SAMPLING_TYPE")))
			fscanf(fptr, "%d", &config.SAMPLING_TYPE);
		else if (!(strcmp(line, "MAX_OBS")))
			fscanf(fptr, "%d", &config.MAX_OBS);

		else if (!(strcmp(line, "EVAL_FLAG")))
			fscanf(fptr, "%d", &config.EVAL_FLAG);
		else if (!(strcmp(line, "EVAL_SEED"))) {
			fscanf(fptr, "%lld", &config.EVAL_SEED[numEvals+1]);
			numEvals++;
			if ( numEvals > maxReps + 1 ) {
				maxReps *= 2;
				config.EVAL_SEED = (long long *) mem_realloc(config.EVAL_SEED, (maxReps+1)*sizeof(long long));
			}
		}
		else if (!(strcmp(line, "EVAL_MIN_ITER")))
			fscanf(fptr, "%d", &config.EVAL_MIN_ITER);
		else if (!(strcmp(line, "EVAL_ERROR")))
			fscanf(fptr, "%lf", &config.EVAL_ERROR);

		else if (!(strcmp(line, "DRO_TYPE"))) {
			int temp;
			fscanf(fptr, "%d", &temp);
			switch (temp) {
			case 0: config.DRO_TYPE = RISK_NEUTRAL; break;
			case 1: config.DRO_TYPE = MOMENT_MATCHING; break;
			default: errMsg("read", "readConfig", "unrecognized parameter in configuration file", 1); break;
			}
		}
		else if (!(strcmp(line, "DRO_PARAM")))
			fscanf(fptr, "%lf", &config.DRO_PARAM);

		else if (!strcmp(line, "//"))
			fgets(comment, 2*BLOCKSIZE, fptr);
		else {
			printf ("%s\n", line);
			errMsg("read", "readConfig", "unrecognized parameter in configuration file", 1);
		}
	}

	fclose(fptr);

	if ( config.MULTIPLE_REP > minimum(numReps, numEvals) ) {
		config.MULTIPLE_REP = minimum(numReps, numEvals);
		printf("Warning: Performing only %d replications due to lack of seeds to run or evaluate.\n", config.MULTIPLE_REP);
	}

	return 0;
}//END readConfig()

int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell, dVector *meanSol) {
	dVector	lb = NULL;
	int 	t;

	/* setup mean value problem which will act as reference for all future computations */
	(*meanSol) = meanProblem(orig, stoc);
	if ( (*meanSol) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to setup and solve mean value problem", 0);
		goto TERMINATE;
	}

	/* decompose the problem into master and subproblem */
	(*prob) = newProb(orig, stoc, tim, lb, config.TOLERANCE);
	if ( (*prob) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to update probType with elements specific to algorithm", 0);
		goto TERMINATE;
	}

	/* ensure that we have a linear programs at all stages */
	t = 0;
	while ( t < tim->numStages ) {
		if ( (*prob)[t++]->sp->type  != PROB_LP )
			printf("Warning :: Stage-%d problem is a mixed-integer program. Solving its linear relaxation.\n", t);
	}

	/* create the cells which will be used in the algorithms */
	(*cell) = newCell(stoc, (*prob), (*meanSol));
	if ( (*cell) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to create the necessary cell structure", 0);
		goto TERMINATE;
	}

	mem_free(lb);
	return 0;
	TERMINATE:
	mem_free(lb);
	return 1;
}//END setupAlgo()

/* This function is used to create cells used in the algorithm */
cellType *newCell(stocType *stoc, probType **prob, dVector xk) {
	cellType    *cell;

	/* Allocate memory to all cells used in the algorithm. */
	if (!(cell = (cellType *) mem_malloc(sizeof(cellType))) )
		errMsg("Memory allocation", "new_cell", "failed to allocate memory to cell",0);
	cell->master = cell->subprob = NULL;
	cell->cuts = NULL;
	cell->omega = NULL; cell->piM = NULL;

	cell->candidX = cell->incumbX = NULL;

	/* stochastic elements */
	cell->omega  = newOmega(stoc);

	/* setup the master problem */
	cell->master = newMaster(prob[0]->sp, prob[0]->lb, cell->omega);
	if ( cell->master == NULL ) {
		errMsg("setup", "newCell", "failed to setup the master problem", 0);
		return NULL;
	}

	/* setup the subproblem */
	cell->subprob = newSubproblem(prob[1]->sp);

	/* Setup the distribution separation problem */
	cell->sep = NULL;

	/* -+-+-+-+-+-+-+-+-+-+-+ Allocating memory to other variables that belongs to cell +-+-+-+-+-+-+-+-+-+- */
	cell->k 	= 0;
	cell->LPcnt = 0;

	/* candidate solution and estimates */
	cell->candidX 	= duplicVector(xk, prob[0]->num->cols);
	cell->candidEst	= prob[0]->lb + vXvSparse(cell->candidX, prob[0]->dBar);

	/* optimality and feasibility flags */
	cell->optFlag 	 = false;
	cell->spFeasFlag = true;
	cell->infeasIncumb = false;

	/* incumbent solution and estimates */
	if (config.MASTER_TYPE == PROB_QP || config.SAMPLING_TYPE == 2 ) {
		cell->incumbX   = duplicVector(xk, prob[0]->num->cols);

		cell->quadScalar= config.MIN_QUAD_SCALAR;     						/* The quadratic scalar, 'sigma'*/
		cell->incumbEst = cell->candidEst;
		cell->incumbChg = false;
		cell->iCutIdx   = 0;
		cell->gamma		= 0.0;

		cell->maxCuts = config.CUT_MULT*prob[0]->num->cols + 1;

		if ( !(cell->piM = (dVector) arr_alloc(prob[0]->num->rows + cell->maxCuts + 1, double)) )
			errMsg("allocation", "newMaster", "cell->piM", 0);
	}
	else {
		cell->incumbX   = NULL;

		cell->quadScalar= 0.0;
		cell->incumbChg = false;
		cell->iCutIdx   = 0;
		cell->gamma		= 0.0;

		cell->maxCuts = config.MAX_ITER;
		cell->piM = NULL;
	}

	cell->cuts = newCuts(cell->maxCuts);

	if ( !(cell->time = (runTime *) mem_malloc(sizeof(runTime)) ) )
		errMsg("setup", "newCell", "cell->runTime", 0);
	cell->time->repTime = 0.0;
	cell->time->iterTime = cell->time->iterAccumTime = 0.0;
	cell->time->masterIter = cell->time->masterAccumTime = 0.0;
	cell->time->subprobIter = cell->time->subprobAccumTime = 0.0;
	cell->time->optTestIter = cell->time->optTestAccumTime = 0.0;
	cell->time->argmaxIter = cell->time->argmaxAccumTime = 0.0;
	cell->time->reduceTime = 0.0;

	/* construct the QP using the current incumbent */
	if ( config.MASTER_TYPE == PROB_QP ) {
		if ( constructQP(prob[0], cell, cell->incumbX, cell->quadScalar) ) {
			errMsg("setup", "newCell", "failed to change the right-hand side after incumbent change", 0);
			return NULL;
		}
	}

	if ( config.SAMPLING_TYPE == 2 ) {
		/* Initialize all the elements which will be used to store dual information */
		cell->lambda = newLambda(cell->omega->cnt*config.MAX_ITER);
		cell->sigma  = newSigma(cell->omega->cnt*config.MAX_ITER);
		cell->delta  = newDelta(cell->omega->cnt*config.MAX_ITER);
	}
	else {
		cell->lambda = NULL;
		cell->sigma  = NULL;
		cell->delta  = NULL;
	}

#if defined(SETUP_CHECK)
	if ( writeProblem(cell->master->lp, "newMaster.lp") ) {
		errMsg("write problem", "new_master", "failed to write master problem to file",0);
		return NULL;
	}
#endif

	return cell;
}//END newCell()

int cleanCellType(cellType *cell, probType *prob, dVector xk) {
	int cnt;

	/* constants and arrays */
	cell->k = 0;
	cell->LPcnt = 0;
	cell->optFlag = false;

	copyVector(xk, cell->candidX, prob->num->cols, true);
	cell->candidEst	= prob->lb + vXvSparse(cell->candidX, prob->dBar);

	if (config.MASTER_TYPE == PROB_QP) {
		copyVector(xk, cell->incumbX, prob->num->cols, true);
		cell->incumbEst = cell->candidEst;
		cell->quadScalar= config.MIN_QUAD_SCALAR;

		cell->iCutIdx = 0;

		cell->incumbChg = true;
	}

	/* oneProblem structures and solver elements */
	for ( cnt = prob->num->rows+cell->cuts->cnt-1; cnt >= prob->num->rows; cnt-- )
		if (  removeRow(cell->master->lp, cnt, cnt) ) {
			errMsg("solver", "cleanCellType", "failed to remove a row from master problem", 0);
			return 1;
		}
	cell->master->mar = prob->num->rows;
	if( changeQPproximal(cell->master->lp, prob->num->cols, cell->quadScalar, 0)) {
		errMsg("algorithm", "cleanCellType", "failed to change the proximal term", 0);
		return 1;
	}

	/* cuts */
	if (cell->cuts) freeCutsType(cell->cuts, true);

	/* stochastic components */
	if ( config.SAMPLING_TYPE == 1 ) {
		cnt = cell->omega->cnt;
		if (cell->omega) freeOmegaType(cell->omega, true);
		cell->omega->cnt = cnt;
	}

	/* reset all the clocks */
	cell->time->repTime = 0.0;
	cell->time->iterTime = cell->time->iterAccumTime = 0.0;
	cell->time->masterIter = cell->time->masterAccumTime = 0.0;
	cell->time->subprobIter = cell->time->subprobAccumTime = 0.0;
	cell->time->optTestIter = cell->time->optTestAccumTime = 0.0;
	cell->time->argmaxIter = cell->time->argmaxAccumTime = 0.0;
	cell->time->reduceTime = 0.0;

	if ( config.MASTER_TYPE == PROB_QP ) {
		if ( constructQP(prob, cell, cell->incumbX, cell->quadScalar) ) {
			errMsg("setup", "newCell", "failed to change the right-hand side after incumbent change", 0);
			return 1;
		}

		cell->incumbChg = false;

#if defined(SETUP_CHECK)
		if ( writeProblem(cell->master->lp, "cleanedQPMaster.lp") ) {
			errMsg("write problem", "new_master", "failed to write master problem to file",0);
			return 1;
		}
#endif
	}

	freeOneProblem(cell->sep);

	return 0;
}//END cleanCellType()

void freeCellType(cellType *cell) {

	if ( cell ) {
		if (cell->master) freeOneProblem(cell->master);
		if (cell->sep) freeOneProblem(cell->sep);
		if (cell->candidX) mem_free(cell->candidX);
		if (cell->incumbX) mem_free(cell->incumbX);
		if (cell->cuts) freeCutsType(cell->cuts, false);
		if (cell->piM) mem_free(cell->piM);
		if (cell->time) mem_free(cell->time);
		if (cell->delta) freeDeltaType(cell->delta, cell->lambda->cnt, cell->omega->cnt);
		if (cell->omega) freeOmegaType(cell->omega, false);
		if (cell->lambda) freeLambdaType(cell->lambda);
		if (cell->sigma) freeSigmaType(cell->sigma);
		mem_free(cell);
	}

}//END freeCellType()
