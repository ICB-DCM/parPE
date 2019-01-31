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
std::unique_ptr<amici::Model> getModel();

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

        file_id = parpe::hdf5CreateFile(outFileArgument.c_str(), true);
        logParPEVersion(file_id);

        // setup data and problem
        dataProvider = std::make_unique<parpe::MultiConditionDataProviderHDF5>(
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

    virtual ~MyOptimizationApplication() override {
        parpe::logProcessStats();
    }

private:
    /** DataProvider as interface to HDF5 data */
    std::unique_ptr<parpe::MultiConditionDataProviderHDF5> dataProvider;
};

int main(int argc, char **argv) {
    MyOptimizationApplication app;
    return app.run(argc, argv);
}
