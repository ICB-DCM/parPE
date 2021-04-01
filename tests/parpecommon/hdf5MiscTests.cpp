#include <gtest/gtest.h>

#include <parpecommon/hdf5Misc.h>

#include "testingMisc.h"

#include <cstdio>
#include <H5Cpp.h>

class HDF5 : public ::testing::Test {

protected:
    void SetUp() override {
        // avoid memory problems
        H5::H5Library::dontAtExit();
        parpe::initHDF5Mutex();

        file = parpe::hdf5CreateFile(tempFileName, true);
    }

    void TearDown() override {
        file.close();
        std::remove(tempFileName.c_str());
    }

    std::string tempFileName {"parpeTest_hdf5Misc.h5"};
    H5::H5File file;
};



TEST_F(HDF5, OpenExistingFileNoOverwrite) {
    EXPECT_THROW(parpe::hdf5CreateFile(tempFileName, false),
                 parpe::HDF5Exception);
}


TEST_F(HDF5, OpenExistingFileOverwrite) {
    file.close();
    file = parpe::hdf5CreateFile(tempFileName, true);
}


TEST_F(HDF5, MutexGetLock) {
    parpe::hdf5MutexGetLock();
}


TEST_F(HDF5, ErrorStackWalker) {
    H5_SAVE_ERROR_HANDLER;

    // provoke error by asking to truncate a file that is already open
    hid_t fileId = H5Fcreate(tempFileName.c_str(), H5F_ACC_TRUNC,
                             H5P_DEFAULT, H5P_DEFAULT);
    EXPECT_TRUE(fileId <= 0);

    auto s = parpe::captureStreamToString([](){
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD,
                 parpe::hdf5ErrorStackWalker_cb, nullptr);
    }, stdout);

    H5_RESTORE_ERROR_HANDLER;

    EXPECT_TRUE(200 < s.size());
}


TEST_F(HDF5, CreateGroup) {
    const char *groupName = "/test";

    EXPECT_FALSE(parpe::hdf5GroupExists(file, groupName));

    parpe::hdf5CreateGroup(file, groupName, false);

    EXPECT_TRUE(parpe::hdf5GroupExists(file, groupName));
}


TEST_F(HDF5, CreateExistingGroup) {
    parpe::hdf5CreateGroup(file, "/test", false);

    H5_SAVE_ERROR_HANDLER;
    EXPECT_THROW(parpe::hdf5CreateGroup(file, "/test", false),
                 parpe::HDF5Exception);
    H5_RESTORE_ERROR_HANDLER;
}


TEST_F(HDF5, EnsureGroupExists) {
    const char *groupName = "/test";

    EXPECT_FALSE(parpe::hdf5GroupExists(file, groupName));

    parpe::hdf5EnsureGroupExists(file, groupName);

    EXPECT_TRUE(parpe::hdf5GroupExists(file, groupName));

    parpe::hdf5EnsureGroupExists(file, groupName);
}

TEST_F(HDF5, StringAttribute) {
    const char *groupName = "/";
    const char *attrName = "testA";
    const char *expAttrValue = "adsf";

    EXPECT_FALSE(parpe::hdf5AttributeExists(file, groupName, attrName));

    parpe::hdf5WriteStringAttribute(file, groupName, attrName, expAttrValue);

    EXPECT_TRUE(parpe::hdf5AttributeExists(file, groupName, attrName));

    H5T_class_t type_class;
    size_t size = 0;
    int ret = H5LTget_attribute_info(file.getId(), groupName, attrName, nullptr,
                                     &type_class, &size);
    EXPECT_TRUE(ret >= 0);
    char actValue[size];

    H5LTget_attribute_string(file.getId(), groupName, attrName, actValue);
    EXPECT_EQ(std::string(expAttrValue), std::string(actValue));
}


TEST_F(HDF5, DatasetDimensions) {
    const char datasetName[] = "bla";
    const int rank = 3;
    const hsize_t dims[rank] = {1,2,3};
    const int buffer[6] = {1};

    EXPECT_FALSE(file.nameExists(datasetName));

    EXPECT_THROW(parpe::hdf5GetDatasetDimensions(
                     file, datasetName, rank,
                     nullptr, nullptr, nullptr, nullptr),
                 H5::Exception);

    EXPECT_TRUE(H5LTmake_dataset_int(file.getId(), datasetName, rank, dims, buffer) >= 0);

    EXPECT_TRUE(file.nameExists(datasetName));

    int d0 = 0, d1 = 0, d2 = 0, d3 = 0;
    parpe::hdf5GetDatasetDimensions(file, datasetName, rank,
                                    &d0, &d1, &d2, &d3);
    EXPECT_EQ((signed)dims[0], d0);
    EXPECT_EQ((signed)dims[1], d1);
    EXPECT_EQ((signed)dims[2], d2);
    EXPECT_EQ(0, d3);
}
