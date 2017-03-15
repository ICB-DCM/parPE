#ifndef OPTIMIZATION_PROBLEM_H
#define OPTIMIZATION_PROBLEM_H


typedef void* (*objectiveFunctionFp)(void *);
typedef void* (*objectiveFunctionGradientFp)(void *);
typedef void* (*intermediateFunctionFp)(void *);

typedef struct OptimizationProblem_tag {
    objectiveFunctionFp objectiveFunction;
    objectiveFunctionGradientFp objectiveFunctionGradient;
    intermediateFunctionFp intermediateFunction;
    double *parametersMin;
    double *parametersMax;
    void *userData;
    char *logFile;
} OptimizationProblem;


#endif
