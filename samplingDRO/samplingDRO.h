/*
 * samplingDRO.h
 *
 *  Created on: Mar 22, 2019
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#ifndef SAMPLINGDRO_H_
#define SAMPLINGDRO_H_

#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"

#define SETUP_CHECK
#define ALGO_CHECK
#undef CUT_CHECK
#undef EVAL_CHECK
#define SEP_CHECK

enum droType {
	RISK_NEUTRAL,
	MOMENT_MATCHING
}droType;

typedef struct{
	int 	MULTIPLE_REP;		/* Number of replications (integer) */
	long long *RUN_SEED;		/* seed used during optimization */

	double 	TOLERANCE; 			/* for zero identity test */

	int		MIN_ITER;			/* minimum number of iterations */
	int		MAX_ITER;			/* maximum number of iterations */
	int		MASTER_TYPE;		/* type of master problem */
	double	EPSILON;			/* Optimality gap */

	int 	TAU;				/* Number of iterations after which incumbent cut is reevaluated */
	int		CUT_MULT;			/* Determines the number of cuts to be used for approximate */

	double	MIN_QUAD_SCALAR;	/* Minimum value for regularizing parameter */
	double 	MAX_QUAD_SCALAR;	/* Maximum value for regularizing parameter */
	double  R1;
	double  R2;
	double  R3;

	int		EVAL_FLAG;
	long long *EVAL_SEED;
	int		EVAL_MIN_ITER;
	double  EVAL_ERROR;

	int		SAA; 				/* Use SAA when continuous distribution in stoch file (1), or not (0) */
	int		SAMPLING_TYPE;		/* Sampling type 0 (full), 1 (SAA), and 2 (Sequential sampling) */
	int		MAX_OBS;			/* Maximum number of iterations before which SAA is invoked */

	int		DRO_TYPE;
	double	DRO_PARAM;
}configType;

typedef struct {
	int		ck;					/* Iteration when the cut was generated */
	int 	numSamples;			/* Number of observations used in computing the cut */
	double  alpha;              /* scalar value for the right-hand side */
	dVector  beta;              /* coefficients of the master problems's primal variables */
	bool	isIncumb;			/* indicates if the cut is an incumbent cut */
	double 	alphaIncumb;		/* right-hand side when using QP master, this is useful for quick updates */
	int 	rowNum;				/* row number for master problem in solver */
	int		omegaID;			/* the observation ID used when multi-cut option is used */
	iVector iStar;				/* Holds the ID for the sigma which is associated with each observation */
	cString	name;
}oneCut;

typedef struct {
	int    	cnt;                /* number of cuts */
	oneCut  **vals;				/* values which define the set of cuts */
}cutsType;

typedef struct {
	double	repTime;
	double 	iterTime;
	double 	masterIter;
	double 	subprobIter;
	double 	optTestIter;
	double  argmaxIter;
	double 	iterAccumTime;
	double 	masterAccumTime;
	double 	subprobAccumTime;
	double 	optTestAccumTime;
	double  argmaxAccumTime;
	double	reduceTime;
}runTime;

typedef struct {
	int		numRV;					/* Number of random variables */
	int 	numObs;					/* Total number of observations */
	int 	cnt;					/* Number of unique observations */
	iVector weights;				/* Frequency of observations */
	dVector	probs;					/* Probability of observation */
	dVector	*vals;					/* Observation values */
	dVector	sampleMean;				/* Sample mean */
} omegaType;

typedef struct {
	int		cnt;					/* number of elements in the structure */
	dVector	*vals;					/* value of duals with random elements in right-hand side */
}lambdaType;

typedef struct{
	double 	pib;					/* scalar pi x b */
	dVector 	piC;					/* dVector pi x C */
} pixbCType;

typedef struct {
	int 		cnt;				/* Number of elements */
	pixbCType 	*vals;				/* product terms */
	iVector		lambdaIdx;			/* Corresponding index in lambdaType */
	iVector		ck;					/* Iteration when the element of generated */
} sigmaType;

typedef struct {
	pixbCType 	**vals;				/* matrix of product terms (rows - entries in lambdaType, columns - entries in omegaType */
} deltaType;

typedef struct {
	int         k;                  /* number of iterations */
	int 		LPcnt; 				/* the number of LPs solved. */

    oneProblem  *master;            /* store master information */
	oneProblem 	*subprob;			/* store subproblem information */

	oneProblem  *sep;				/* distribution separation problem */

	dVector     candidX;            /* primal solution of the master problem */
	double      candidEst;          /* objective value master problem */

	dVector     incumbX;			/* incumbent master solution */
	double      incumbEst;			/* estimate at incumbent solution */
	double		gamma;				/* Improvement in objective function */
	double 		quadScalar; 		/* the proximal parameter/quadratic scalar 'sigma' */
	bool        incumbChg;			/* set to be true if the incumbent solution has changed in an iteration */
	int     	iCutIdx;			/* index of incumbent cuts in cell->cuts structure. */
	int         iCutUpdt;			/* iteration when incumbent cut was last updated. */
	dVector		piM;				/* Dual vector for the master problem (original and the cuts) */

    int      	maxCuts;            /* maximum number of cuts to be used*/
	cutsType    *cuts;              /* optimality cuts */

	bool        optFlag;			/* Flag to indicate optimality of a cell */
	bool		spFeasFlag;			/* Subproblem feasibility flag */
	bool        infeasIncumb;		/* Indicates if the incumbent solution is infeasible */

	runTime		*time;				/* Run time structure */

	omegaType 	*omega;				/* all realizations observed during the algorithm */
	lambdaType	*lambda;			/* Duals corresponding to rows with randomness */
	sigmaType	*sigma;				/* Deterministic components of cut coefficients */
	deltaType	*delta;				/* Stochastic components of cut coefficients */
}cellType;

/* samplingDRO.c */
void parseCmdLine(int argc, char *argv[], cString *probName, cString *inputDir);
void printHelpMenu();
void writeOptimizationStatistics(FILE *soln, FILE *incumb, probType **prob, cellType *cell, int rep);
void printOptimizationSummary(cellType *cell);

/* setup.c */
int readConfig();
int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell, dVector *meanSol);
cellType *newCell(stocType *stoc, probType **prob, dVector xk);
int cleanCellType(cellType *cell, probType *prob, dVector xk);
void freeCellType(cellType *cell);

/* algo.c */
int algo(oneProblem *orig, timeType *tim, stocType *stoc, cString inputDir, cString probName);
int solveFixedDROCell(probType **prob, cellType *cell);
int solveDRSDCell(stocType *stoc, probType **prob, cellType *cell);

/* optimal.c */
bool optimal(probType *prob, cellType *cell);

/* master.c */
int solveMaster(numType *num, sparseVector *dBar, cellType *cell);
int checkImprovement(probType *prob, cellType *cell, int candidCut);
int replaceIncumbent(cellType *cell, double candidEst, int numCols);
int constructQP(probType *prob, cellType *cell, dVector incumbX, double quadScalar);
int changeEtaCol(LPptr lp, int numRows, int numCols, int k, cutsType *cuts);
int changeQPproximal(LPptr lp, int numCols, double sigma, int numEta);
int changeQPrhs(probType *prob, cellType *cell, dVector xk);
int changeQPbds(LPptr lp, int numCols, dVector bdl, dVector bdu, dVector xk);
int addCut2Master(cellType *cell, cutsType *cuts, oneCut *cut, int lenX, int obsID);
oneProblem *newMaster(oneProblem *orig, double lb, omegaType *omega);

/* subproblem.c */
int solveSubprob(probType *prob, oneProblem *subproblem, dVector Xvect, dVector obsVals, bool *spFeasFlag, double *subprobTime, dVector piS, double *mubBar);
dVector computeRHS(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, dVector X, dVector obs);
dVector computeCostCoeff(numType *num, coordType *coord, sparseVector *dBar, dVector obs);
int computeMU(LPptr lp, int numCols, double *mubBar);
oneProblem *newSubproblem(oneProblem *subprob);
void chgRHSwSoln(sparseVector *bBar, sparseMatrix *Cbar, dVector rhs, dVector X);
int chgRHSwObserv(LPptr lp, numType *num, coordType *coord, dVector observ, dVector spRHS, dVector X);
int chgObjxwObserv(LPptr lp, numType *num, coordType *coord, dVector cost, iVector indices, dVector observ);

/* cuts.c */
int formDeterministicCut(probType *prob, cellType *cell, dVector Xvect);
int formStochasticCut(probType *prob, cellType *cell, dVector Xvect, int obsStar, bool newOmegaFlag);
oneCut *newCut(int numX, int currentIter, int numObs);
cutsType *newCuts(int maxCuts);
double maxCutHeight(cutsType *cuts, int currIter, dVector xk, int betaLen, double lb, bool scale);
double cutHeight(oneCut *cut, int currIter, dVector xk, int betaLen, double lb, bool scale);
int computeIstar(numType *num, coordType *coord, sigmaType *sigma, deltaType *delta, dVector piCbarX, dVector Xvect,
		int obs, double *argmax);
int reduceCuts();
void freeOneCut(oneCut *cut);
void freeCutsType(cutsType *cuts, bool partial);

/* stocUpdate.c */
int stochasticUpdates(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, lambdaType *lambda, sigmaType *sigma,
		deltaType *delta, omegaType *omega, int iter, dVector pi, double mubBar);
lambdaType *newLambda(int numLambda);
sigmaType *newSigma(int numSigma);
deltaType *newDelta(int numDelta);
int calcOmega(dVector observ, dVector stocMean, omegaType *omega, bool *newOmegaFlag, double TOLERANCE);
void refineOmega(omegaType *omega, dVector stocMean);
int calcLambda(numType *num, coordType *coord, dVector Pi, lambdaType *lambda, bool *newLambdaFlag);
int calcSigma(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *CBar, dVector pi, double mubBar,
		int idxLambda, bool newLambdaFlag, int currentIter, sigmaType *sigma);
int calcDeltaRow(numType *num, coordType *coord, omegaType *omega, lambdaType *lambda, int lambdaIdx, deltaType *delta);
void freeLambdaType (lambdaType *lambda);
void freeSigmaType (sigmaType *sigma);
void freeDeltaType (deltaType *delta, int numLambda, int numOmega);
omegaType *newOmega(stocType *stoc);
void freeOmegaType(omegaType *omega, bool partial);

/* separation.c */
int obtainProbDist(oneProblem *sep, dVector stocMean, omegaType *omega, dVector spObj, int obsStar, bool newOmegaFlag);
int updtDistSepProb_MM(oneProblem *sep, dVector stocMean, omegaType *omega, dVector spObj, int obsStar, bool newOmegaFlag, int numMoments);
oneProblem *newDistSepProb(dVector stocMean, omegaType *omega);
oneProblem *newDistSepProb_MM(dVector stocMean, omegaType *omega, int numMoments);

/* evalution.c */
int evaluate(FILE *soln, stocType *stoc, probType **prob, cellType *cell, dVector Xvect);
void printEvaluationSummary(FILE *soln, double mean, double stdev, int cnt);
void writeEvaluationStatistics(FILE *soln, double mean, double stdev, int cnt);

#endif /* SAMPLINGDRO_H_ */
