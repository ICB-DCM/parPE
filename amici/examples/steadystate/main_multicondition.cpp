#include "steadyStateMultiConditionDataprovider.h"
#include "wrapfunctions.h"

#include <optimizationOptions.h>
#include <loadBalancerWorker.h>
#include <hdf5Misc.h>
#include <logging.h>
#include <optimizationApplication.h>
#include <misc.h>
#include <hdf5Misc.h>

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
 * @brief The SteadystateApplication class subclasses parpe::OptimizationApplication
 * which provides a frame for a standalone program to solve a multi-start local optimization
 * problem.
 */

class SteadystateApplication : public parpe::OptimizationApplication {
  public:
    using OptimizationApplication::OptimizationApplication;

    ~SteadystateApplication() override = default;

    void initProblem(std::string inFileArgument,
                             std::string outFileArgument) override {

        // The same file should only be opened/created once, an then only be reopened
        file_id = parpe::hdf5CreateFile(outFileArgument.c_str(), true);
        logParPEVersion(file_id);

        dataProvider = std::make_unique<SteadyStateMultiConditionDataProvider>(
                    getModel(), inFileArgument);

        // read options from file
        auto optimizationOptions = parpe::OptimizationOptions::fromHDF5(dataProvider->getHdf5FileId());

        // Create one instance for the problem, one for the application for clear ownership
        auto multiCondProb = new parpe::MultiConditionProblem(
                    dataProvider.get(), &loadBalancer,
                    std::make_unique<parpe::Logger>(),
                    // TODO remove this resultwriter
                    std::make_unique<parpe::OptimizationResultWriter>(
                        file_id,
                        std::string("/multistarts/"))
                    );

        // hierarchical optimization?
        if(optimizationOptions->hierarchicalOptimization) {
            problem.reset(new parpe::HierarchicalOptimizationProblemWrapper(
                              std::unique_ptr<parpe::MultiConditionProblem>(multiCondProb),
                              dataProvider.get())
                          );
        } else {
            problem.reset(multiCondProb);
        }

        problem->setOptimizationOptions(*optimizationOptions);

        // On master, copy input data to result file
        if(parpe::getMpiRank() < 1)
            dataProvider->copyInputData(file_id);

        auto ms = new parpe::MultiConditionProblemMultiStartOptimizationProblem(
                    dataProvider.get(),
                    problem->getOptimizationOptions(),
                    multiCondProb->getResultWriter(),
                    &loadBalancer,
                    std::make_unique<parpe::Logger>()
                    );
        multiStartOptimizationProblem.reset(ms);
    }

    std::unique_ptr<SteadyStateMultiConditionDataProvider> dataProvider;
};

/**
 * @brief The SteadystateLocalOptimizationApplication class overrides the multi-start optimization
 * in the base class and performs only a single optimization run. This is mostly for debugging.
 */
class SteadystateLocalOptimizationApplication : public SteadystateApplication {
  public:

    using SteadystateApplication::SteadystateApplication;

    void runMaster() override {
        // Single optimization
        getLocalOptimum(problem.get());
    }
};


int main(int argc, char **argv) {
    int status = EXIT_SUCCESS;

    // SteadystateLocalOptimizationApplication app(argc, argv);
    SteadystateApplication app;
    status = app.run(argc, argv);

    return status;
}


