#ifndef PARAMETERMAPPER_H
#define PARAMETERMAPPER_H

class ParameterMapper {
  public:
    ParameterMapper();

    ParameterMapper(numCommonParameters, numConditionSpecificParamsPerCondition,
                    numConditions);

    int getSimulationParameterIndex(int optimizationParameterIndex);

    int getOptimizationParameterIndex(int simulationParameterIndex);

  protected:
    int numCommonParameters = 0;
    int numConditionSpecificParamsPerCondition = 0;
    int numConditions = 0;
};

#endif // PARAMETERMAPPER_H
