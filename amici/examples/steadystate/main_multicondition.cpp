#include "steadyStateMultiConditionDataprovider.h"
#include "optimizationOptions.h"
#include "wrapfunctions.h"
#include <LoadBalancerWorker.h>
#include <hdf5Misc.h>
#include <logging.h>
#include <multiConditionProblemResultWriter.h>
#include <optimizationApplication.h>
#include <misc.h>
#include <iostream>
#include <unistd.h>

/** @file
 *
 * This example demonstrates the use of the loadbalancer / queue for
 * parallel ODE simulation.
 * The example is based on the `steadystate` example included in AMICI.
 *
 * This model has 5 parameters, 3 states, 4 condition specific fixed parameters.
 *
 * To run, e.g.: mpiexec -np 4
 * ../parPE-build/amici/examples/steadystate/example_steadystate_multi -o
 * steadystate_`date +%F` amici/examples/steadystate/data.h5
 */

/**
 * @brief The SteadyStateMultiConditionProblem class subclasses parpe::MultiConditionProblem
 * to set some problem specific options like initial parameters
 */
class SteadyStateMultiConditionProblem : public parpe::MultiConditionProblem {
  public:
    SteadyStateMultiConditionProblem(
        SteadyStateMultiConditionDataProvider *dp, parpe::LoadBalancerMaster *loadBalancer)
        : MultiConditionProblem(dp, loadBalancer) {

        std::unique_ptr<parpe::OptimizationOptions> options(parpe::OptimizationOptions::fromHDF5(
                                                                 dataProvider->getHdf5FileId()));

        optimizationOptions = *options.get();
        std::fill(initialParameters_.begin(), initialParameters_.end(), 0);
    }
};

/**
 * @brief The SteadystateApplication class subclasses parpe::OptimizationApplication
 * which provides a frame for a standalone program to solve a multi-start local optimization
 * problem.
 */

class SteadystateApplication : public parpe::OptimizationApplication {
  public:
    using OptimizationApplication::OptimizationApplication;

    virtual void initProblem(std::string inFileArgument,
                             std::string outFileArgument) override {
        model = std::unique_ptr<Model>(getModel());
        dataProvider = std::make_unique<SteadyStateMultiConditionDataProvider>(model.get(), inFileArgument);

        problem = std::unique_ptr<parpe::MultiConditionProblem>(
                    new SteadyStateMultiConditionProblem(dataProvider.get(), &loadBalancer));

        parpe::JobIdentifier id;
        resultWriter = std::make_unique<parpe::MultiConditionProblemResultWriter>(outFileArgument, true, id);

        problem->resultWriter = std::make_unique<parpe::MultiConditionProblemResultWriter>(*resultWriter);
        problem->resultWriter->setJobId(id);

        if(parpe::getMpiRank() < 1)
            dataProvider->copyInputData(resultWriter->file_id);

    }

    std::unique_ptr<SteadyStateMultiConditionDataProvider> dataProvider;
    std::unique_ptr<Model> model;
};

/**
 * @brief The SteadystateLocalOptimizationApplication class overrides the multi-start optimization
 * in the base class and performs only a single optimization run. This is mostly for debugging.
 */
class SteadystateLocalOptimizationApplication : public SteadystateApplication {
  public:

    using SteadystateApplication::SteadystateApplication;

    virtual int runMaster() override {
        // Single optimization
        return getLocalOptimum(problem.get());
    }
};


int main(int argc, char **argv) {
    int status = EXIT_SUCCESS;

    // SteadystateLocalOptimizationApplication app(argc, argv);
    SteadystateApplication app;
    status = app.run(argc, argv);

    return status;
}


