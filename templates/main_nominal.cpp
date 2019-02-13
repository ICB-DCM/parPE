#include <parpeamici/multiConditionDataProvider.h>
#include <parpeamici/multiConditionProblem.h>
#include <parpecommon/logging.h>
#include <parpecommon/misc.h>
#include <parpeoptimization/optimizationOptions.h>
#include <amici/model.h>
#include <iostream>
#include <memory>

// to avoid including model-specific header files
std::unique_ptr<amici::Model> getModel();
using namespace parpe;

int main(int argc, char **argv) {
#ifndef NDEBUG
    // Set stdout to unbuffered when debugging
    setbuf(stdout, NULL);
#endif

    std::string inFileArgument = "";

    // Check command line arguments
    if (argc != 2) {
        fprintf(stderr, "Error: must provide HDF5 input file as first and only "
                        "argument.\n");
        return 1;
    } else {
        inFileArgument = argv[1];
    }


    parpe::logmessage(parpe::LOGLVL_INFO,
                      "Reading options and data from '%s'.",
                      inFileArgument.c_str());

    // setup data and problem
    MultiConditionDataProviderHDF5 dataProvider(getModel(), inFileArgument);
    auto options = OptimizationOptions::fromHDF5(dataProvider.getHdf5FileId());

    MultiConditionProblem problem {&dataProvider};

    // We assume we have nominal parameter stored as starting point 0
    // TODO: should be stored separately
    int start_idx = 0;
    auto optimizationParams = options->getStartingPoint(dataProvider.getHdf5FileId(), start_idx);

    double fval = NAN;
    std::vector<double> gradient(optimizationParams.size(), NAN);
    problem.costFun->evaluate(optimizationParams, fval, gradient);

    std::cout<<gradient<<std::endl;
    std::for_each(gradient.begin(), gradient.end(), [](double &d){ d = std::fabs(d); });
    auto maxAbsGradient = *std::max_element(gradient.begin(), gradient.end());
    std::for_each(gradient.begin(), gradient.end(), [fval](double &d){ d /= std::fabs(fval); });
    auto maxRelGradient = *std::max_element(gradient.begin(), gradient.end());
    parpe::logmessage(LOGLVL_INFO, "Max(abs(grad)) = " + std::to_string(maxAbsGradient));
    parpe::logmessage(LOGLVL_INFO, "Max(abs(grad)/fval) = " + std::to_string(maxRelGradient));

    parpe::logmessage(LOGLVL_INFO, "Likelihood: " + std::to_string(fval));
}
