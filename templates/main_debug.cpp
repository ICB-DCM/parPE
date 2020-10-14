#include <parpeamici/multiConditionDataProvider.h>
#include <parpecommon/logging.h>
#include <parpecommon/misc.h>
#include <parpeoptimization/optimizationOptions.h>
#include <amici/model.h>
#include <iostream>
#include <memory>

// to avoid including model-specific header files
namespace amici::generic_model {
std::unique_ptr<amici::Model> getModel();
}
using namespace parpe;

int main(int argc, char **argv) {
#ifndef NDEBUG
    // Set stdout to unbuffered when debugging
    setbuf(stdout, NULL);
#endif

    std::string inFileArgument = "/home/dweindl/src/benchmarkProblem/20190205221009_Speedy_v4_Jan2019_generic_degradation_r415549/Speedy_v4_Jan2019_generic_degradation_r415549.bak.h5";

    parpe::logmessage(parpe::LOGLVL_INFO,
                      "Reading options and data from '%s'.",
                      inFileArgument.c_str());

    // setup data and problem
    MultiConditionDataProviderHDF5 dataProvider(
        amici::generic_model::getModel(), inFileArgument);
    auto options = OptimizationOptions::fromHDF5(dataProvider.getHdf5FileId());


    auto model = dataProvider.getModel();
    //model->setTimepoints({1e3});
    model->requireSensitivitiesForAllParameters();

    int condition_idx = 0;
    int start_idx = 0;

    auto optimizationParams = options->getStartingPoint(dataProvider.getHdf5FileId(), start_idx);
    dataProvider.updateSimulationParametersAndScale(condition_idx, optimizationParams, *model);
    auto edata = dataProvider.getExperimentalDataForCondition(condition_idx);

    auto solver = dataProvider.getSolver();
    solver->setSensitivityOrder(amici::SensitivityOrder::none);
    solver->setSensitivityMethod(amici::SensitivityMethod::adjoint);

    solver->setSensitivityOrder(amici::SensitivityOrder::first);
    solver->setMaxSteps(10000);
    solver->setMaxStepsBackwardProblem(10000);

    WallTimer timer;
    auto rdata = amici::runAmiciSimulation(*solver, edata.get(), *model);
    std::cout<<timer.getTotal()<<" s"<<std::endl;

    std::cout<<"Numsteps: "<<rdata->numsteps<<std::endl;
    std::cout<<"NumstepsB: "<<rdata->numstepsB<<std::endl;
    //logProcessStats();
}
