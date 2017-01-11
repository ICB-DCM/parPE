#ifndef OBJECTIVE_FUNCTION_H
#define OBJECTIVE_FUNCTION_H

#define OBJECTIVE_FUNCTION_H_DEBUG_LEVEL 0

#include<include/udata.h>
#include<include/edata.h>
#include <include/rdata.h>
#include<dataprovider.h>

/**
 * @brief evaluateObjectiveFunction Evaluate objective function and gradient (optional)
 * @param theta Parameter vector
 * @param lenTheta Length of parameter vector
 * @param path Experimental data to use
 * @param objectiveFunctionValue Output: pointer to buffer for objective function value
 * @param objectiveFunctionGradient Output: pointer to buffer for objective function gradient or NULL to not evaluate gradient
 * @param scaling Parameter scaling
 * @return 0 if successful
 */

int evaluateObjectiveFunction(const double theta[], int lenTheta, datapath path, double *objectiveFunctionValue, double *objectiveFunctionGradient, AMI_parameter_scaling scaling);

/**
 * @brief getSteadystateSolution Simulate the model until steady state
 * @param udata: Model data and options
 * @param edata: ExpData TODO: could remove in case of no events?
 * @param status: return status, 0 if successful
 * @param iterationDone: Iterations of simulation to t=10^9 until steady-state is reached
 * @return
 */

ReturnData *getSteadystateSolution(UserData *udata, ExpData *edata, int *status, int *iterationDone);

/**
 * @brief getSteadystateSolutionForExperiment same as evaluateObjectiveFunction, but edata is obtained automatically and optionally returned
 * @param path
 * @param udata
 * @param status
 * @param _edata Output: the requested ExpData
 * @param iterationsDone
 * @return
 */

ReturnData *getSteadystateSolutionForExperiment(datapath path, UserData *udata, int *status, ExpData **_edata, int *iterationsDone);

#endif
