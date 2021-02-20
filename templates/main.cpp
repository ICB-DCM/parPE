#include <parpeamici/optimizationApplication.h>
#include <parpeamici/multiConditionDataProvider.h>
#include <parpeamici/multiConditionProblem.h>
#include <parpeoptimization/optimizationResultWriter.h>
#include <parpecommon/logging.h>
#include <parpeoptimization/optimizationOptions.h>
#include <parpecommon/misc.h>

#include <amici/model.h>

#include <memory>

// to avoid including model-specific header files
namespace amici::generic_model {
    std::unique_ptr<amici::Model> getModel();
}

class MyOptimizationApplication : public parpe::OptimizationApplication {
public:
    using OptimizationApplication::OptimizationApplication;

    virtual void initProblem(std::string inFileArgument,
                             std::string outFileArgument) override
    {
        if (!isWorker())
            parpe::logmessage(parpe::LOGLVL_INFO,
                              "Reading options and data from '%s'.",
                              inFileArgument.c_str());

        auto h5Outfile = parpe::hdf5CreateFile(outFileArgument, true);
        logParPEVersion(h5Outfile);

        // setup data and problem
        dataProvider = std::make_unique<parpe::MultiConditionDataProviderHDF5>(
            amici::generic_model::getModel(), inFileArgument);

        // read options from file
        auto h5Infile = dataProvider->getHdf5File();
        auto optimizationOptions = parpe::OptimizationOptions::fromHDF5(h5Infile);

        // Create one instance for the problem, one for the application for clear ownership
        auto multiCondProb = new parpe::MultiConditionProblem(
                    dataProvider.get(), &loadBalancer,
                    std::make_unique<parpe::Logger>(),
                    // TODO remove this resultwriter
                    std::make_unique<parpe::OptimizationResultWriter>(
                        h5Outfile,
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
            dataProvider->copyInputData(h5Outfile);

        auto ms = new parpe::MultiConditionProblemMultiStartOptimizationProblem(
                    dataProvider.get(),
                    problem->getOptimizationOptions(),
                    multiCondProb->getResultWriter(),
                    &loadBalancer,
                    std::make_unique<parpe::Logger>()
                    );
        multiStartOptimizationProblem.reset(ms);
    }

    virtual ~MyOptimizationApplication() override {
        parpe::logProcessStats();
    }

private:
    /** DataProvider as interface to HDF5 data */
    std::unique_ptr<parpe::MultiConditionDataProviderHDF5> dataProvider;
};

int main(int argc, char **argv) {
#ifndef NDEBUG
    // Set stdout to unbuffered when debugging
    setbuf(stdout, NULL);
#endif
    MyOptimizationApplication app;
    return app.run(argc, argv);
}
