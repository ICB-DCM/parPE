#include "multiConditionProblemResultWriter.h"


MultiConditionProblemResultWriter::MultiConditionProblemResultWriter() : OptimizationResultWriter()
{

}

MultiConditionProblemResultWriter::MultiConditionProblemResultWriter(hid_t file_id, JobIdentifier id) : OptimizationResultWriter(file_id)
{
    this->id = id;
    rootPath = getIterationPath();

}

MultiConditionProblemResultWriter::MultiConditionProblemResultWriter( const char *filename, bool overwrite, JobIdentifier id) : OptimizationResultWriter(filename, overwrite)
{
    this->id = id;
    rootPath = getIterationPath();
}

std::string MultiConditionProblemResultWriter::getIterationPath()
{

    char fullGroupPath[1024];
    sprintf(fullGroupPath, "/multistarts/%d/iteration/%d/",
            id.idxLocalOptimization, id.idxLocalOptimizationIteration);
    return std::string(fullGroupPath);
}

std::string MultiConditionProblemResultWriter::getSimulationPath()
{
    char fullGroupPath[1024];
    sprintf(fullGroupPath, "/multistarts/%d/iteration/%d/genotype/%d/experiment/%d/",
            id.idxLocalOptimization, id.idxLocalOptimizationIteration, 0, id.idxConditions);
    return std::string(fullGroupPath);
}

void MultiConditionProblemResultWriter::logSimulation(const double *theta, double llh, const double *gradient, double timeElapsedInSeconds, int nTheta, int numStates, double *states, double *stateSensi, double *y, int jobId, int iterationsUntilSteadystate)
{
    std::string pathStr = getSimulationPath();
    const char *fullGroupPath = pathStr.c_str();


    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "negLogLikelihood", &llh, 1);
    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath, "jobId", &jobId, 1);

    if(gradient)
        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "negLogLikelihoodGradient", gradient, nTheta);

    if(theta)
        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "p", theta, nTheta);

    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "evalFTime", &timeElapsedInSeconds, 1);

    if(states)
        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "X", states, numStates);

    if(y)
        hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath, "Y", y, 1);

    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath, "iterationsUntilSteadystate", &iterationsUntilSteadystate, 1);

    if(stateSensi)
        hdf5CreateOrExtendAndWriteToDouble3DArray(file_id, fullGroupPath, "sX", stateSensi, numStates, nTheta);

    flushResultWriter();

}
