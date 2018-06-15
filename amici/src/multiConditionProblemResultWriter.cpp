#include "multiConditionProblemResultWriter.h"
#include <logging.h>
#include <cmath>

namespace parpe {

MultiConditionProblemResultWriter::MultiConditionProblemResultWriter()
    : OptimizationResultWriter() {}

MultiConditionProblemResultWriter::MultiConditionProblemResultWriter(
    hid_t file_id)
    : OptimizationResultWriter(file_id)
{
    auto lock = hdf5MutexGetLock();
    this->file_id = H5Freopen(file_id);
}

MultiConditionProblemResultWriter::MultiConditionProblemResultWriter(
    hid_t file_id, JobIdentifier id)
    : MultiConditionProblemResultWriter(file_id)
{
    setJobId(id);
}

MultiConditionProblemResultWriter::MultiConditionProblemResultWriter(
    std::string filename, bool overwrite, JobIdentifier id)
    : OptimizationResultWriter(filename, overwrite)
{
    auto lock = hdf5MutexGetLock();
    setJobId(id);
}

MultiConditionProblemResultWriter::MultiConditionProblemResultWriter(const MultiConditionProblemResultWriter &other)
    : MultiConditionProblemResultWriter(other.file_id, other.id)
{

}

void MultiConditionProblemResultWriter::logLocalOptimizerObjectiveFunctionEvaluation(
        gsl::span<const double> parameters, double objectiveFunctionValue,
        gsl::span<const double>objectiveFunctionGradient,
        int numIterations, int numFunctionCalls, double timeElapsedInSeconds)
{
    OptimizationResultWriter::logLocalOptimizerObjectiveFunctionEvaluation(
                parameters, objectiveFunctionValue,
                objectiveFunctionGradient, numIterations, numFunctionCalls, timeElapsedInSeconds);

    if(std::isinf(objectiveFunctionValue) || std::isnan(objectiveFunctionValue))
        printObjectiveFunctionFailureMessage();

}

std::string MultiConditionProblemResultWriter::getIterationPath(int iterationIdx) const {

    char fullGroupPath[1024];
    sprintf(fullGroupPath, "/multistarts/%d/iteration/%d/",
            id.idxLocalOptimization, iterationIdx);
    return fullGroupPath;
}


void MultiConditionProblemResultWriter::printObjectiveFunctionFailureMessage() const {
    char strBuf[100];
    id.sprint(strBuf);
    logmessage(LOGLVL_ERROR, "%s: Objective function evaluation failed!",
               strBuf);
}

void MultiConditionProblemResultWriter::saveLocalOptimizerResults(double finalNegLogLikelihood,
                                                   gsl::span<const double> optimalParameters, double masterTime,
                                                   int exitStatus) const
{
    char strBuf[100];
    id.sprint(strBuf);
    logmessage(LOGLVL_INFO, "%s: Optimizer status %d, final llh: %e, time: %f.",
               strBuf, exitStatus, finalNegLogLikelihood, masterTime);

    OptimizationResultWriter::saveLocalOptimizerResults(finalNegLogLikelihood, optimalParameters,
                              masterTime, exitStatus);
}

std::string MultiConditionProblemResultWriter::getOptimizationPath() const {
    char fullGroupPath[1024];
    sprintf(fullGroupPath, "/multistarts/%d/", id.idxLocalOptimization);
    return fullGroupPath;
}

void MultiConditionProblemResultWriter::setJobId(JobIdentifier id) {
    this->id = id;
    setRootPath(getIterationPath(id.idxLocalOptimizationIteration));
}

JobIdentifier MultiConditionProblemResultWriter::getJobId() { return id; }

} // namespace parpe
