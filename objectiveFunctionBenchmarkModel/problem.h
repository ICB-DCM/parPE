#ifndef PROBLEM_H
#define PROBLEM_H

#include "optimizationProblem.h"
#include "resultwriter.h"
#include "dataprovider.h"
#include "include/amici.h"

typedef struct {
    int nTheta;
    double *theta;
    UserData udata;
    ExpData edata;
    ReturnData rdata;
    double *gradient;
    double objectiveFunctionValue;
    Datapath datapath;
    int scaling;
} MyUserData;


OptimizationProblem *getBenchmarkOptimizationProblem(Datapath path);

void freeBenchmarkProblem(OptimizationProblem *problem);
#endif
