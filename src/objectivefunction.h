#ifndef OBJECTIVE_FUNCTION_H
#define OBJECTIVE_FUNCTION_H

#define OBJECTIVE_FUNCTION_H_DEBUG_LEVEL 0

#include<include/amici.h>
#include<dataprovider.h>
#include<masterqueue.h>


// data to be returned to queuer
typedef struct queueResults_tag {
    int datasetID;
    int lenllh;
    double *llh;
    int lensllh;
    double *sllh;
    int lenx;
    double *x;
    int leny;
    double *y;
} queueResults;

queueData *getDataForSimulation(datapath datapath);

int evaluateObjectiveFunction(const double theta[], int lenTheta, datapath path, double *objectiveFunctionValue, double *objectiveFunctionGradient, AMI_parameter_scaling scaling);

ReturnData *getSteadystateSolution(UserData *udata, ExpData *edata, int *status, int *iterationDone);

ReturnData *getSteadystateSolutionForExperiment(datapath path, UserData *udata, int *status, ExpData **_edata, int *iterationsDone);


#endif
