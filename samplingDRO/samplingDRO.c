/*
 * samplingDRO.c
 *
 *  Created on: Mar 22, 2019
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "samplingDRO.h"

cString outputDir;
long long	MEM_USED = 0;	/* amount of memory allocated each iteration */
configType config;

int main(int argc, char * argv[]) {
	cString inputDir = NULL, probName = NULL;
	oneProblem *orig = NULL;
	timeType *tim = NULL;
	stocType *stoc = NULL;
	outputDir = NULL;

	/* read algorithm configuration file */
	if ( readConfig() )
		goto TERMINATE;

	/* read problem information */
	parseCmdLine(argc, argv, &probName, &inputDir);

	/* read problem SMPS input files */
	if ( readFiles(inputDir, probName, &orig, &tim, &stoc) ) {
		errMsg("read", "main", "failed to read problem files using SMPS reader", 0);
		goto TERMINATE;
	}

	/* set up output directory: using the outputDir in configuration file and the input problem name */
	createOutputDir(outputDir, "samplingDRO", probName);

	/* launch the algorithm */
	if ( algo(orig, tim, stoc, inputDir, probName) ) {
		errMsg("allocation", "main", "failed to solve the problem using 2-SD algorithm", 0);
		goto TERMINATE;
	}

	TERMINATE:
	mem_free(inputDir); mem_free(probName); mem_free(outputDir);
	return 0;
}//END main()

void parseCmdLine(int argc, char *argv[], cString *probName, cString *inputDir) {

	for(int i=1; (i < argc); i++) {
		if ( argv[i][0] == '-' ) {
			switch ((argv[i])[1]) {
			case '?': printHelpMenu(); exit(0);
			case 'p': {
				(*probName) = (cString) arr_alloc(2*BLOCKSIZE, char);
				strcpy((*probName), argv[++i]); break;
			}
			case 'i': {
				(*inputDir) = (cString) arr_alloc(2*BLOCKSIZE, char);
				strcpy((*inputDir), argv[++i]); break;
			}
			case 'o': {
				outputDir = (cString) arr_alloc(2*BLOCKSIZE, char);
				strcpy(outputDir, argv[++i]); break;
			}
		}}
		else {
			printf("Input options must begin with a '-'. Use '-?' for help.\n"); exit(0);
		}
	}

	if ( probName == NULL || inputDir == NULL || outputDir == NULL ) {
		printf("Problem name, input and output directory are mandatory input.\n");
		if ( (*probName) ) mem_free((*probName));
		if ( outputDir ) mem_free(outputDir);
		if ( (*inputDir) ) mem_free((*inputDir));
		closeSolver(); exit(0);
	}

}//END parseCmdLine()

/* We allow only a few of the parameters to be selected through command line input. */
void printHelpMenu() {

	printf("Input options:\n");
	/* Problem name, input and output directory */
	printf("         -p string  -> problem name.\n");
	printf("         -i string  -> input directory where the problem SMPS files are saved.\n");
	printf("         -o string  -> output directory where the result files will be written.\n");

}//END helpMenu()
