#include "multiConditionProblemResultWriter.h"


MultiConditionProblemResultWriter::MultiConditionProblemResultWriter() : OptimizationResultWriter()
{

}

MultiConditionProblemResultWriter::MultiConditionProblemResultWriter(hid_t file_id, JobIdentifier id) : OptimizationResultWriter(file_id)
{
    setJobId(id);
}

MultiConditionProblemResultWriter::MultiConditionProblemResultWriter( const char *filename, bool overwrite, JobIdentifier id) : OptimizationResultWriter(filename, overwrite)
{
    setJobId(id);
}

std::string MultiConditionProblemResultWriter::getIterationPath(int iterationIdx)
{

    char fullGroupPath[1024];
    sprintf(fullGroupPath, "/multistarts/%d/iteration/%d/",
            id.idxLocalOptimization, id.idxLocalOptimizationIteration);
    return fullGroupPath;
}

std::string MultiConditionProblemResultWriter::getSimulationPath(JobIdentifier id)
{
    char fullGroupPath[1024];
    sprintf(fullGroupPath, "/multistarts/%d/iteration/%d/condition/%d/",
            id.idxLocalOptimization, id.idxLocalOptimizationIteration, id.idxConditions);
    return fullGroupPath;
}

std::string MultiConditionProblemResultWriter::getOptimizationPath()
{
    char fullGroupPath[1024];
    sprintf(fullGroupPath, "/multistarts/%d/", id.idxLocalOptimization);
    return fullGroupPath;
}

void MultiConditionProblemResultWriter::setJobId(JobIdentifier id)
{
    this->id = id;
    rootPath = getIterationPath(id.idxLocalOptimizationIteration);
}

JobIdentifier MultiConditionProblemResultWriter::getJobId()
{
    return id;
}

void MultiConditionProblemResultWriter::logSimulation(JobIdentifier id, const double *theta, double llh, const double *gradient, double timeElapsedInSeconds, int nTheta, int numStates, double *states, double *stateSensi, double *y, int jobId, int iterationsUntilSteadystate)
{
    std::string pathStr = getSimulationPath(id);
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

void MultiConditionProblemResultWriter::saveLocalOptimizerResults(double finalNegLogLikelihood, const double *optimalParameters, int numParameters, double masterTime, int exitStatus)
{
    std::string oldPath = rootPath;
    rootPath = getOptimizationPath();
    OptimizationResultWriter::saveLocalOptimizerResults(finalNegLogLikelihood, optimalParameters, numParameters, masterTime, exitStatus);
    rootPath = oldPath;
}
