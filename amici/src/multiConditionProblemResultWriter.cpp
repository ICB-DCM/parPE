#include "multiConditionProblemResultWriter.h"
#include <cmath>

namespace parpe {

MultiConditionProblemResultWriter::MultiConditionProblemResultWriter()
    : OptimizationResultWriter() {}

MultiConditionProblemResultWriter::MultiConditionProblemResultWriter(
    hid_t file_id)
    : OptimizationResultWriter(file_id) {
}

MultiConditionProblemResultWriter::MultiConditionProblemResultWriter(
    hid_t file_id, JobIdentifier id)
    : OptimizationResultWriter(file_id) {
    setJobId(id);
}

MultiConditionProblemResultWriter::MultiConditionProblemResultWriter(
    std::string filename, bool overwrite, JobIdentifier id)
    : OptimizationResultWriter(filename, overwrite) {
    setJobId(id);
}

MultiConditionProblemResultWriter::MultiConditionProblemResultWriter(const MultiConditionProblemResultWriter &other) : MultiConditionProblemResultWriter(other.file_id, other.id)
{

}

std::string MultiConditionProblemResultWriter::getIterationPath(int iterationIdx) {

    char fullGroupPath[1024];
    sprintf(fullGroupPath, "/multistarts/%d/iteration/%d/",
            id.idxLocalOptimization, id.idxLocalOptimizationIteration);
    return fullGroupPath;
}

std::string
MultiConditionProblemResultWriter::getSimulationPath(JobIdentifier id) {
    char fullGroupPath[1024];
    sprintf(fullGroupPath, "/multistarts/%d/iteration/%d/condition/%d/",
            id.idxLocalOptimization, id.idxLocalOptimizationIteration,
            id.idxConditions);
    return fullGroupPath;
}

std::string MultiConditionProblemResultWriter::getOptimizationPath() {
    char fullGroupPath[1024];
    sprintf(fullGroupPath, "/multistarts/%d/", id.idxLocalOptimization);
    return fullGroupPath;
}

void MultiConditionProblemResultWriter::setJobId(JobIdentifier id) {
    this->id = id;
    setRootPath(getIterationPath(id.idxLocalOptimizationIteration));
}

JobIdentifier MultiConditionProblemResultWriter::getJobId() { return id; }

void MultiConditionProblemResultWriter::logSimulation(
    JobIdentifier id, const double *theta, double llh, const double *gradient,
    double timeElapsedInSeconds, int nTheta, int numStates, double *states,
    double *stateSensi, int numY, double *y, int jobId,
    int iterationsUntilSteadystate, int status) {
    // TODO replace by SimulationResultWriter
    if (!gradient && !logLineSearch)
        return;

    std::string pathStr = getSimulationPath(id);
    const char *fullGroupPath = pathStr.c_str();

    hdf5CreateOrExtendAndWriteToDouble2DArray(
        file_id, fullGroupPath, "simulationLogLikelihood", &llh, 1);

    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath, "jobId",
                                           &jobId, 1);

    if (gradient) {
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file_id, fullGroupPath, "simulationLogLikelihoodGradient", gradient,
            nTheta);
    } else {
        double dummyGradient[nTheta];
        std::fill_n(dummyGradient, nTheta, NAN);
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file_id, fullGroupPath, "simulationLogLikelihoodGradient",
            dummyGradient, nTheta);
    }

    if (theta)
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file_id, fullGroupPath, "simulationParameters", theta, nTheta);

    hdf5CreateOrExtendAndWriteToDouble2DArray(file_id, fullGroupPath,
                                              "simulationWallTimeInSec",
                                              &timeElapsedInSeconds, 1);

    if (states)
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file_id, fullGroupPath, "simulationStates", states, numStates);

    if (y)
        hdf5CreateOrExtendAndWriteToDouble2DArray(
            file_id, fullGroupPath, "simulationObservables", y, numY);

    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath,
                                           "iterationsUntilSteadystate",
                                           &iterationsUntilSteadystate, 1);

    if (stateSensi)
        hdf5CreateOrExtendAndWriteToDouble3DArray(
            file_id, fullGroupPath, "simulationStateSensitivities", stateSensi,
            numStates, nTheta);

    hdf5CreateOrExtendAndWriteToInt2DArray(file_id, fullGroupPath,
                                           "simulationStatus", &status, 1);

    flushResultWriter();
}

} // namespace parpe
