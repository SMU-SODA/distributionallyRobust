# Decomposition Algorithms for Distributionally Robust Optimization

This repository consists of a suite of decomposition algorithms, primarily based on the notion of sequential sampling. The sequential sampling algorithms find their genesis in the stochastic decomposition algorithm of Higle and Sen (1991). The repository currently consists of the following algorithms. 
1. Distributionally robust stochastic decomposition:

    Details can be found in Gangammanavar, H and Bansal, M (2020). Stochastic Decomposition Method for Two-Stage Distributionally Robust Optimization, available on Optimization Online [(here)](http://www.optimization-online.org/DB_HTML/2020/11/8091.html)

The code is developed by Harsha Gangammanavar [(profile](https://github.com/gjharsha)/[website)](https://gjharsha.github.io/).

## Input file format
The implementation of accepts stochastic programs described using the SMPS file format:

* Core file: Describes the optimization problem corresponding to a single scenario
* Time file: Provides information regarding decomposition of the problem into master and subproblem
* Stoch file: Describes the stochastic information associated with the problem. (Currently, INDEP and BLOCK types are supported)

For more information about the SMPS file format [here](https://doi.org/10.1137/1.9780898718799.ch2)

## Support
Please report bugs [here on GitHub](https://github.com/SMU-SODA/distributionallyRobust/issues).

## Installation
Note: Only Unix systems have been tested.
### Prerequisite: 
  * CPLEX should be available on your machine. 
  * Git, if you want to directly use the commands below.
  * This repository uses certain utilities and data sets from the `spAlgorithms` respository. 

### Steps
  1. Download the SD source codes.  
    * `git clone https://github.com/SMU-SODA/distributionallyRobust.git`  
  2. Compile the source files using appropriate tools that suit your system.   
  3. Setup a directory to write output files.  
    * `mkdir spOutput`  
  5. Execute the algorithm. The input options for the algorithm are:  
         Input options:  
             `-p` string  -> problem name.  
             `-i` string  -> input directory where the problem SMPS files are saved.  
             `-o` string  -> output directory where the result files will be written.  
     * Example: `-p pgp2 -i ../../spAlgorithms/spInput/ -o ../../spOutput
