#include <gtest/gtest.h>

#include <parpecommon/hdf5Misc.h>

#include "testingMisc.h"

#include <cstdio>
#include <H5Cpp.h>

class hdf5Misc : public ::testing::Test {

protected:
    void SetUp() override {
        // avoid memory problems
        H5::H5Library::dontAtExit();
        parpe::initHDF5Mutex();

        fileId = parpe::hdf5CreateFile(std::tmpnam(tempFileName), false);
    }

    void TearDown() override {
        if(fileId)
            H5Fclose(fileId);
        std::remove(tempFileName);
    }

    char tempFileName[TMP_MAX];
    hid_t fileId = 0;
};



TEST_F(hdf5Misc, testOpenExistingFileNoOverwrite) {
    EXPECT_THROW(parpe::hdf5CreateFile(tempFileName, false),
                 parpe::HDF5Exception);
}


TEST_F(hdf5Misc, testOpenExistingFileOverwrite) {
    H5Fclose(fileId);
    fileId = parpe::hdf5CreateFile(tempFileName, true);
}


TEST_F(hdf5Misc, testMutexGetLock) {
    parpe::hdf5MutexGetLock();
}


TEST_F(hdf5Misc, testErrorStackWalker) {
    H5_SAVE_ERROR_HANDLER;

    // provoke error by asking to truncate a file that is already open
    hid_t fileId = H5Fcreate(tempFileName, H5F_ACC_TRUNC,
                             H5P_DEFAULT, H5P_DEFAULT);
    EXPECT_TRUE(fileId <= 0);

    auto s = parpe::captureStreamToString([](){
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD,
                 parpe::hdf5ErrorStackWalker_cb, nullptr);
    }, stdout);

    H5_RESTORE_ERROR_HANDLER;

    EXPECT_TRUE(200 < s.size());
}


TEST_F(hdf5Misc, testCreateGroup) {
    const char *groupName = "/test";

    EXPECT_FALSE(parpe::hdf5GroupExists(fileId, groupName));

    parpe::hdf5CreateGroup(fileId, groupName, false);

    EXPECT_TRUE(parpe::hdf5GroupExists(fileId, groupName));
}


TEST_F(hdf5Misc, testCreateExistingGroup) {
    parpe::hdf5CreateGroup(fileId, "/test", false);

    H5_SAVE_ERROR_HANDLER;
    EXPECT_THROW(parpe::hdf5CreateGroup(fileId, "/test", false),
                 parpe::HDF5Exception);
    H5_RESTORE_ERROR_HANDLER;
}


TEST_F(hdf5Misc, testEnsureGroupExists) {
    const char *groupName = "/test";

    EXPECT_FALSE(parpe::hdf5GroupExists(fileId, groupName));

    parpe::hdf5EnsureGroupExists(fileId, groupName);

    EXPECT_TRUE(parpe::hdf5GroupExists(fileId, groupName));

    parpe::hdf5EnsureGroupExists(fileId, groupName);
}

TEST_F(hdf5Misc, testStringAttribute) {
    const char *groupName = "/";
    const char *attrName = "testA";
    const char *expAttrValue = "adsf";

    EXPECT_FALSE(parpe::hdf5AttributeExists(fileId, groupName, attrName));

    parpe::hdf5WriteStringAttribute(fileId, groupName, attrName, expAttrValue);

    EXPECT_TRUE(parpe::hdf5AttributeExists(fileId, groupName, attrName));

    H5T_class_t type_class;
    size_t size = 0;
    int ret = H5LTget_attribute_info(fileId, groupName, attrName, nullptr,
                                     &type_class, &size);
    EXPECT_TRUE(ret >= 0);
    char actValue[size];

    H5LTget_attribute_string(fileId, groupName, attrName, actValue);
    EXPECT_EQ(std::string(expAttrValue), std::string(actValue));
}


TEST_F(hdf5Misc, testDatasetDimensions) {
    const char datasetName[] = "bla";
    const int rank = 3;
    const hsize_t dims[rank] = {1,2,3};
    const int buffer[6] = {1};

    EXPECT_FALSE(parpe::hdf5DatasetExists(fileId, datasetName));

    EXPECT_THROW(parpe::hdf5GetDatasetDimensions(
                     fileId, datasetName, rank,
                     nullptr, nullptr, nullptr, nullptr),
                 H5::Exception);

    EXPECT_TRUE(H5LTmake_dataset_int(fileId, datasetName, rank, dims, buffer) >= 0);

    EXPECT_TRUE(parpe::hdf5DatasetExists(fileId, datasetName));

    int d0 = 0, d1 = 0, d2 = 0, d3 = 0;
    parpe::hdf5GetDatasetDimensions(fileId, datasetName, rank,
                                    &d0, &d1, &d2, &d3);
    EXPECT_EQ((signed)dims[0], d0);
    EXPECT_EQ((signed)dims[1], d1);
    EXPECT_EQ((signed)dims[2], d2);
    EXPECT_EQ(0, d3);
}

// TODO:
// hdf5CreateExtendableDouble2DArray
// hdf5CreateOrExtendAndWriteToDouble2DArray
// hdf5CreateOrExtendAndWriteToDouble3DArray
// hdf5CreateOrExtendAndWriteToInt2DArray
// hdf5CreateExtendableInt2DArray
// hdf5CreateExtendableDouble3DArray
// hdf5Extend2ndDimensionAndWriteToDouble2DArray
// hdf5Extend3rdDimensionAndWriteToDouble3DArray
// hdf5Extend2ndDimensionAndWriteToInt2DArray
// hdf5Read2DDoubleHyperslab
// hdf5Read3DDoubleHyperslab
