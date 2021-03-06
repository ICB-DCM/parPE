#include <gtest/gtest.h>

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

// need prototype here, otherwise mess with headers
// (including ceres.h causes some errors with EIGEN)
namespace parpe {
void setCeresOption(const std::pair<const std::string, const std::string> &pair,
                    ceres::GradientProblemSolver::Options* options);
} // namespace parpe
#endif


TEST(OptimizationOptions, setGetOptionStr) {
    parpe::OptimizationOptions o;
    std::string key = "str";
    std::string expVal = "testStr";
    o.setOption(key, expVal);
    auto actVal = o.getStringOption(key);

    EXPECT_EQ(expVal, actVal);
}

TEST(OptimizationOptions, setGetOptionInt) {
    parpe::OptimizationOptions o;
    std::string key = "str";
    double expVal = 1.23;
    o.setOption(key, expVal);
    auto actVal = o.getDoubleOption(key);

    EXPECT_NEAR(expVal, actVal, 1e-15);
}

TEST(OptimizationOptions, setGetOptionDouble) {
    parpe::OptimizationOptions o;
    std::string key = "str";
    auto expVal = 123;
    o.setOption(key, expVal);
    auto actVal = o.getIntOption(key);

    EXPECT_EQ(expVal, actVal);
}


TEST(OptimizationOptions, getNonExistingOption) {
    parpe::OptimizationOptions o;

    EXPECT_THROW(o.getIntOption("missingKey"), std::invalid_argument);
}

#ifdef PARPE_ENABLE_IPOPT
TEST(OptimizationOptions, setIpOptOptions) {
    std::string key = "max_iter";
    int expVal = 10;

    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    auto options = app->Options();

    parpe::OptimizationOptions o;
    o.setOption(key, expVal);
    o.for_each<Ipopt::SmartPtr<Ipopt::OptionsList>*>(parpe::setIpOptOption,
                                                     &options);

    int actVal = 0;
    EXPECT_EQ(true, options->GetIntegerValue(key, actVal, ""));
    EXPECT_EQ(expVal, actVal);
}
#endif

#ifdef PARPE_ENABLE_CERES
TEST(OptimizationOptions, setCeresOptions) {
    std::string key = "max_num_iterations";
    int expVal = 10;

    ceres::GradientProblemSolver::Options options;

    parpe::OptimizationOptions o;
    o.setOption(key, expVal);
    o.for_each<ceres::GradientProblemSolver::Options*>(parpe::setCeresOption,
                                                       &options);

    int actVal = options.max_num_iterations;

    // NOTE: disabled for ceres 1.11 compatibility
    // CHECK_TRUE(options.IsValid(nullptr));
    EXPECT_EQ(expVal, actVal);
}
#endif

TEST(OptimizationOptions, fromHDF5) {
    const char* tmpName = "parpeTest_fromHDF5.h5";
    auto _ = gsl::finally([tmpName] { remove(tmpName); });

    // fail on non-existing file (hide hdf5 errors)
    parpe::captureStreamToString([tmpName](){
        EXPECT_THROW(parpe::OptimizationOptions::fromHDF5(tmpName),
                     parpe::HDF5Exception);
    }, stdout);

    // create file
    auto file = parpe::hdf5CreateFile(tmpName, false);
    parpe::hdf5CreateGroup(file, "/optimizationOptions/ceres", true);
    int optimizer = 1;
    H5LTset_attribute_int(file.getId(), "/optimizationOptions", "optimizer",
                          &optimizer, 1);
    H5LTset_attribute_int(file.getId(), "/optimizationOptions/ceres", "someOption",
                          &optimizer, 1);
    hsize_t dims[] = {2, 3};
    double buf[] = {1, 2, 3,
                    4, 5, 6};
    H5LTmake_dataset_double(file.getId(), "/optimizationOptions/randomStarts", 2,
                            dims, buf);

    auto startingPoint = parpe::OptimizationOptions::getStartingPoint(file, 0);
    EXPECT_EQ(1, startingPoint[0]);
    EXPECT_EQ(4, startingPoint[1]);

    auto o = parpe::OptimizationOptions::fromHDF5(tmpName);
    EXPECT_EQ(optimizer, static_cast<int>(o->optimizer));
    EXPECT_EQ(optimizer, o->getIntOption("someOption"));
    EXPECT_TRUE(o->toString().size() > 50);
}
