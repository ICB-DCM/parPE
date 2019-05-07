#include <parpeamici/simulationResultWriter.h>

#include <parpecommon/hdf5Misc.h>

#include <cstdio>
#include <cmath> // NAN
#include <utility>

namespace parpe {

SimulationResultWriter::SimulationResultWriter(H5::H5File const& file,
                                               std::string rootPath)
    : rootPath(std::move(rootPath))
{
    auto lock = hdf5MutexGetLock();
    this->file = file;

    updatePaths();
}


SimulationResultWriter::SimulationResultWriter(const std::string &hdf5FileName,
                                               std::string rootPath)
    : rootPath(std::move(rootPath))
{
    auto lock = hdf5MutexGetLock();

    H5_SAVE_ERROR_HANDLER;

    // Open for append or create
    try {
        file = H5::H5File(hdf5FileName, H5F_ACC_RDWR);
    } catch (H5::FileIException const&) {
        // create if doesn't exist
        file = H5::H5File(hdf5FileName, H5F_ACC_EXCL);
    }
    H5_RESTORE_ERROR_HANDLER;

    updatePaths();
}


void SimulationResultWriter::createDatasets(hsize_t numSimulations)
{

    auto lock = parpe::hdf5MutexGetLock();


    double fillValueDbl = NAN;   /* Fill value for the dataset */
    H5::DSetCreatPropList propList;
    propList.setFillValue(H5::PredType::NATIVE_DOUBLE, &fillValueDbl);

    hdf5EnsureGroupExists(file, rootPath);

    if(saveX || saveYMes || saveYSim)
        hdf5EnsureGroupExists(file, timePath);
    if(saveX)
        hdf5EnsureGroupExists(file, xPath);
    if(saveYMes)
        hdf5EnsureGroupExists(file, yMesPath);
    if(saveYSim)
        hdf5EnsureGroupExists(file, ySimPath);

    // Individual datasets will be created per condition, since we need
    // condition-specific number of timepoints
    // Can only create llh dataset in advance
    if(saveLlh && !hdf5DatasetExists(file.getId(), llhPath.c_str())) {
        hsize_t dims[] = {numSimulations};
        H5::DataSpace dataspace(1, dims);
        hsize_t one = 1;
        propList.setChunk(1, &one);
        file.createDataSet(llhPath, H5::PredType::NATIVE_DOUBLE,
                           dataspace, propList);
    }

    file.flush(H5F_SCOPE_LOCAL);
}


void SimulationResultWriter::createDatasets(int numberOfSimulations)
{
    createDatasets(static_cast<hsize_t>(numberOfSimulations));
}


void SimulationResultWriter::saveSimulationResults(
        const amici::ExpData *edata,
        const amici::ReturnData *rdata,
        int simulationIdx)
{
    saveMeasurements(edata->getObservedData(), edata->nt(),
                     edata->nytrue(), simulationIdx);
    saveTimepoints(rdata->ts, simulationIdx);
    saveModelOutputs(rdata->y,  edata->nt(), edata->nytrue(), simulationIdx);
    saveStates(rdata->x, rdata->nt, rdata->nx, simulationIdx);
    saveLikelihood(rdata->llh, simulationIdx);

    auto lock = parpe::hdf5MutexGetLock();

    file.flush(H5F_SCOPE_LOCAL);
}

void SimulationResultWriter::saveTimepoints(gsl::span<const double> timepoints,
                                            int simulationIdx) const
{
    if(timepoints.empty())
        return;

    auto lock = parpe::hdf5MutexGetLock();

    // Create dataset
    constexpr int rank = 1;
    hsize_t dims[rank] = {static_cast<hsize_t>(timepoints.size())};
    H5::DataSpace dataspace(rank, dims);
    auto dataset = file.createDataSet(
                timePath + "/" + std::to_string(simulationIdx),
                H5::PredType::NATIVE_DOUBLE, dataspace);

    dataset.write(timepoints.data(), H5::PredType::NATIVE_DOUBLE);
}


void SimulationResultWriter::saveMeasurements(
        gsl::span<const double> measurements, int nt,
        int nytrue, int simulationIdx) const
{
    if(!saveYMes || nt < 1 || nytrue < 1 || measurements.empty()) {
        return;
    }

    RELEASE_ASSERT(measurements.size() ==
                   static_cast<decltype(measurements)::index_type>(nt * nytrue),
                   "");

    auto lock = parpe::hdf5MutexGetLock();

    // Create dataset
    constexpr int rank = 2;
    hsize_t dims[rank] = {static_cast<hsize_t>(nt),
                          static_cast<hsize_t>(nytrue)};
    H5::DataSpace dataspace(rank, dims);
    auto dataset = file.createDataSet(
                yMesPath + "/" + std::to_string(simulationIdx),
                H5::PredType::NATIVE_DOUBLE, dataspace);

    dataset.write(measurements.data(), H5::PredType::NATIVE_DOUBLE);
}

void SimulationResultWriter::saveModelOutputs(
        gsl::span<const double> outputs, int nt,
        int nytrue, int simulationIdx) const
{
    if(!saveYSim || nt < 1 || nytrue < 1 || outputs.empty()) {
        return;
    }

    RELEASE_ASSERT(outputs.size() ==
                   static_cast<decltype(outputs)::index_type>(nt * nytrue), "");

    auto lock = parpe::hdf5MutexGetLock();

    // Create dataset
    constexpr int rank = 2;
    hsize_t dims[rank] = {static_cast<hsize_t>(nt),
                          static_cast<hsize_t>(nytrue)};
    H5::DataSpace dataspace(rank, dims);
    auto dataset = file.createDataSet(
                ySimPath + "/" + std::to_string(simulationIdx),
                H5::PredType::NATIVE_DOUBLE, dataspace);

    dataset.write(outputs.data(), H5::PredType::NATIVE_DOUBLE);
}


void SimulationResultWriter::saveStates(
        gsl::span<const double> states, int nt, int nx, int simulationIdx) const
{
    if(!saveX || states.empty() || nx < 1 || nt < 1) {
        return;
    }

    RELEASE_ASSERT(states.size() ==
                   static_cast<decltype(states)::index_type>(nt * nx), "");

    auto lock = parpe::hdf5MutexGetLock();

    // Create dataset
    constexpr int rank = 2;
    hsize_t dims[rank] = {static_cast<hsize_t>(nt),
                          static_cast<hsize_t>(nx)};
    H5::DataSpace dataspace(rank, dims);
    auto dataset = file.createDataSet(
                xPath + "/" + std::to_string(simulationIdx),
                H5::PredType::NATIVE_DOUBLE, dataspace);

    dataset.write(states.data(), H5::PredType::NATIVE_DOUBLE);
}


void SimulationResultWriter::saveLikelihood(double llh, int simulationIdx) const
{
    if(!saveLlh) {
        return;
    }

    auto lock = parpe::hdf5MutexGetLock();

    auto dataset = file.openDataSet(llhPath);

    auto filespace = dataset.getSpace();
    hsize_t count[] = {1};
    hsize_t start[] = {static_cast<hsize_t>(simulationIdx)};
    filespace.selectHyperslab(H5S_SELECT_SET, count, start);

    auto memspace = dataset.getSpace();
    hsize_t start2[] = {0};
    memspace.selectHyperslab(H5S_SELECT_SET, count, start2);

    dataset.write(&llh, H5::PredType::NATIVE_DOUBLE, memspace, filespace);
}


H5::H5File SimulationResultWriter::reopenFile()
{
    auto lock = hdf5MutexGetLock();
    return H5::H5File(file.getId());
}


void SimulationResultWriter::updatePaths()
{
    yMesPath = rootPath + "/yMes";
    ySimPath = rootPath + "/ySim";
    xPath = rootPath + "/x";
    llhPath  = rootPath + "/llh";
    timePath  = rootPath + "/t";
}


} // namespace parpe
