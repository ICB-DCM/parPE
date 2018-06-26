#include "simulationResultWriter.h"
#include <hdf5Misc.h>

#include <cstdio>
#include <cmath> // NAN

namespace parpe {

SimulationResultWriter::SimulationResultWriter(H5::H5File const& file, std::string const& rootPath)
    : rootPath(rootPath)
{
    auto lock = hdf5MutexGetLock();
    this->file = file;
    updatePaths();
}

SimulationResultWriter::SimulationResultWriter(const std::string &hdf5FileName, std::string const& rootPath)
    : rootPath(rootPath)
{
    auto lock = hdf5MutexGetLock();

    H5_SAVE_ERROR_HANDLER;

    try {
        file = H5::H5File(hdf5FileName, H5F_ACC_RDWR);
    } catch (H5::FileIException const& e) {
        // create if doesn't exist
        file = H5::H5File(hdf5FileName, H5F_ACC_EXCL);
    }
    H5_RESTORE_ERROR_HANDLER;

    updatePaths();
}


void SimulationResultWriter::createDatasets(const amici::Model &model,
                                            int numberOfSimulations)
{
    auto numSimulations = static_cast<hsize_t>(numberOfSimulations);
    auto ny = static_cast<hsize_t>(model.nytrue);
    auto nx = static_cast<hsize_t>(model.nxtrue);
    auto nt = static_cast<hsize_t>(model.nt());

    auto lock = parpe::hdf5MutexGetLock();


    double fillValueDbl = NAN;   /* Fill value for the dataset */
    H5::DSetCreatPropList propList;
    propList.setFillValue(H5::PredType::NATIVE_DOUBLE, &fillValueDbl);

    parpe::hdf5EnsureGroupExists(file.getId(), rootPath.c_str());

    if((saveYMes || saveYSim) && model.nt() > 0 && model.nytrue > 0) {
        // observables
        constexpr int rank = 3;
        hsize_t dims[rank] = {numSimulations, nt, ny};
        H5::DataSpace dataspace(rank, dims);
        hsize_t chunkDims[rank] = {1, dims[1], dims[2]};
        propList.setChunk(rank, chunkDims);

        if(saveYMes && !hdf5DatasetExists(file.getId(), yMesPath.c_str())) // TODO should check dimensions if exists
            file.createDataSet(yMesPath, H5::PredType::NATIVE_DOUBLE, dataspace, propList);
        if(saveYSim && !hdf5DatasetExists(file.getId(), ySimPath.c_str()))
            file.createDataSet(ySimPath, H5::PredType::NATIVE_DOUBLE, dataspace, propList);
    }

    if(saveX && !hdf5DatasetExists(file.getId(), xPath.c_str())) {
        int rank = 3;
        hsize_t dims[] = {numSimulations, nt, nx};
        H5::DataSpace dataspace(rank, dims);
        hsize_t chunkDims[] = {1, dims[1], dims[2]};
        propList.setChunk(rank, chunkDims);
        file.createDataSet(xPath, H5::PredType::NATIVE_DOUBLE, dataspace, propList);
    }

    if(saveLlh && !hdf5DatasetExists(file.getId(), llhPath.c_str())) {
        hsize_t dims[] = {numSimulations};
        H5::DataSpace dataspace(1, dims);
        hsize_t one = 1;
        propList.setChunk(1, &one);
        file.createDataSet(llhPath, H5::PredType::NATIVE_DOUBLE, dataspace, propList);
    }

    file.flush(H5F_SCOPE_LOCAL);
}


void SimulationResultWriter::saveSimulationResults(
        const amici::ExpData *edata,
        const amici::ReturnData *rdata,
        int simulationIdx)
{
    auto simulationIdxH = static_cast<hsize_t>(simulationIdx);
    auto ny = static_cast<hsize_t>(rdata->ny);
    auto nx = static_cast<hsize_t>(rdata->nx);
    auto nt = static_cast<hsize_t>(rdata->nt);

    auto lock = parpe::hdf5MutexGetLock();

    if(saveYMes && edata->nt > 0 && edata->nytrue > 0) {
        auto dataset = file.openDataSet(yMesPath);

        auto filespace = dataset.getSpace();
        hsize_t count[] = {1, nt, ny};
        hsize_t start[] = {simulationIdxH, 0, 0};
        filespace.selectHyperslab(H5S_SELECT_SET, count, start);

        auto memspace = dataset.getSpace();
        hsize_t start2[] = {0, 0, 0};
        memspace.selectHyperslab(H5S_SELECT_SET, count, start2);

        dataset.write(edata->my.data(), H5::PredType::NATIVE_DOUBLE, memspace, filespace);
    }

    if(saveYSim && edata->nt > 0 && edata->nytrue > 0) {
        auto dataset = file.openDataSet(ySimPath);

        auto filespace = dataset.getSpace();
        hsize_t count[] = {1, nt, ny};
        hsize_t start[] = {simulationIdxH, 0, 0};
        filespace.selectHyperslab(H5S_SELECT_SET, count, start);

        auto memspace = dataset.getSpace();
        hsize_t start2[] = {0, 0, 0};
        memspace.selectHyperslab(H5S_SELECT_SET, count, start2);

        dataset.write(rdata->y.data(), H5::PredType::NATIVE_DOUBLE, memspace, filespace);
    }

    if(saveX) {
        auto dataset = file.openDataSet(xPath);

        auto filespace = dataset.getSpace();
        hsize_t count[] = {1, nt, nx};
        hsize_t start[] = {simulationIdxH, 0, 0};
        filespace.selectHyperslab(H5S_SELECT_SET, count, start);

        auto memspace = dataset.getSpace();
        hsize_t start2[] = {0, 0, 0};
        memspace.selectHyperslab(H5S_SELECT_SET, count, start2);

        dataset.write(rdata->x.data(), H5::PredType::NATIVE_DOUBLE, memspace, filespace);
    }

    if(saveLlh) {
        auto dataset = file.openDataSet(llhPath);

        auto filespace = dataset.getSpace();
        hsize_t count[] = {1};
        hsize_t start[] = {simulationIdxH};
        filespace.selectHyperslab(H5S_SELECT_SET, count, start);

        auto memspace = dataset.getSpace();
        hsize_t start2[] = {0};
        memspace.selectHyperslab(H5S_SELECT_SET, count, start2);

        dataset.write(&rdata->llh, H5::PredType::NATIVE_DOUBLE, memspace, filespace);
    }

    file.flush(H5F_SCOPE_LOCAL);
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

}



} // namespace parpe
