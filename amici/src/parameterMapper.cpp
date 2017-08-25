#include "parameterMapper.h"

ParameterMapper::ParameterMapper() {}

ParameterMapper::ParameterMapper(numCommonParameters,
                                 numConditionSpecificParamsPerCondition,
                                 numConditions) {}

int ParameterMapper::getSimulationParameterIndex(
    int optimizationParameterIndex) {
    abort(); //
}

int ParameterMapper::getOptimizationParameterIndex(int simulationParameterIndex,
                                                   int conditionIndex) {
    return numCommonParameters +
           numConditionSpecificParamsPerCondition * conditionIndex;
}
