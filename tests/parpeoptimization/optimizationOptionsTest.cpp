#include <CppUTest/TestHarness.h>
#include <CppUTestExt/MockSupport.h>

#include <iostream>

#include <parpecommon/parpeConfig.h>
#include <parpeoptimization/optimizationOptions.h>
#include <parpecommon/hdf5Misc.h>

#include "../parpecommon/testingMisc.h"


#ifdef PARPE_ENABLE_IPOPT
#include <parpeoptimization/localOptimizationIpopt.h>
#include <parpeoptimization/localOptimizationIpoptTNLP.h>
#endif

#ifdef PARPE_ENABLE_CERES
#include <parpeoptimization/localOptimizationCeres.h>
#include <ceres/gradient_problem_solver.h>

// need prototype here, otherwise mess with headers (including ceres.h causes some errors with EIGEN)
namespace parpe {
void setCeresOption(const std::pair<const std::string, const std::string> &pair, ceres::GradientProblemSolver::Options* options);
} // namespace parpe
#endif

// clang-format off
TEST_GROUP(optimizationOptions){
    void setup() {
        mock().clear();
    }

    void teardown() {
        mock().checkExpectations();
        mock().clear();
    }
};
// clang-format on


TEST(optimizationOptions, setGetOptionStr) {
    parpe::OptimizationOptions o;
    std::string key = "str";
    std::string expVal = "testStr";
    o.setOption(key, expVal);
    auto actVal = o.getStringOption(key);

    CHECK_EQUAL(expVal, actVal);
}

TEST(optimizationOptions, setGetOptionInt) {
    parpe::OptimizationOptions o;
    std::string key = "str";
    double expVal = 1.23;
    o.setOption(key, expVal);
    auto actVal = o.getDoubleOption(key);

    DOUBLES_EQUAL(expVal, actVal, 1e-15);
}

TEST(optimizationOptions, setGetOptionDouble) {
    parpe::OptimizationOptions o;
    std::string key = "str";
    auto expVal = 123;
    o.setOption(key, expVal);
    auto actVal = o.getIntOption(key);

    CHECK_EQUAL(expVal, actVal);
}


TEST(optimizationOptions, getNonExistingOption) {
    parpe::OptimizationOptions o;

    CHECK_THROWS(std::invalid_argument, o.getIntOption("missingKey"));
}

#ifdef PARPE_ENABLE_IPOPT
TEST(optimizationOptions, setIpOptOptions) {
    std::string key = "max_iter";
    int expVal = 10;

    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    auto options = app->Options();

    parpe::OptimizationOptions o;
    o.setOption(key, expVal);
    o.for_each<Ipopt::SmartPtr<Ipopt::OptionsList>*>(parpe::setIpOptOption, &options);

    int actVal = 0;
    CHECK_EQUAL(true, options->GetIntegerValue(key, actVal, ""));
    CHECK_EQUAL(expVal, actVal);
}
#endif

#ifdef PARPE_ENABLE_CERES
TEST(optimizationOptions, setCeresOptions) {
    std::string key = "max_num_iterations";
    int expVal = 10;

    ceres::GradientProblemSolver::Options options;

    parpe::OptimizationOptions o;
    o.setOption(key, expVal);
    o.for_each<ceres::GradientProblemSolver::Options*>(parpe::setCeresOption, &options);

    int actVal = options.max_num_iterations;

    // NOTE: disabled for ceres 1.11 compatibility CHECK_TRUE(options.IsValid(nullptr));
    CHECK_EQUAL(expVal, actVal);
}
#endif

TEST(optimizationOptions, fromHDF5) {
    char tmpName[TMP_MAX];
    if(!std::tmpnam(tmpName))
        std::abort();

    // TODO: hide hdf5 errors
//    parpe::captureStreamToString([tmpName](){
    H5_SAVE_ERROR_HANDLER;
        CHECK_THROWS(parpe::HDF5Exception, parpe::OptimizationOptions::fromHDF5(tmpName));
    H5_RESTORE_ERROR_HANDLER;
//    }, stdout);

    // create file
    hid_t fileId = parpe::hdf5CreateFile(tmpName, false);
    parpe::hdf5CreateGroup(fileId, "/optimizationOptions/ceres", true);
    int optimizer = 1;
    H5LTset_attribute_int(fileId, "/optimizationOptions", "optimizer", &optimizer, 1);
    H5LTset_attribute_int(fileId, "/optimizationOptions/ceres", "someOption", &optimizer, 1);
    hsize_t dims[] = {2, 3};
    double buf[] = {1, 2, 3,
                    4, 5, 6};
    H5LTmake_dataset_double(fileId, "/optimizationOptions/randomStarts", 2, dims, buf);

    // TODO: throw instead of pront
    auto startingPoint = parpe::OptimizationOptions::getStartingPoint(fileId, 0);
    CHECK_EQUAL(1, startingPoint[0]);
    CHECK_EQUAL(4, startingPoint[1]);
    H5Fclose(fileId);

    auto o = parpe::OptimizationOptions::fromHDF5(tmpName);
    CHECK_EQUAL(optimizer, static_cast<int>(o->optimizer));
    CHECK_EQUAL(optimizer, o->getIntOption("someOption"));
    CHECK_TRUE(o->toString().size() > 50);

    remove(tmpName);
}
