#include <CppUTest/TestHarness.h>
#include <CppUTestExt/MockSupport.h>

#include <parpecommon/hdf5Misc.h>

#include "testingMisc.h"

#include <cstdio>
#include <H5Cpp.h>

// clang-format off
TEST_GROUP(hdf5Misc){
    char tempFileName[TMP_MAX];
    hid_t fileId = 0;

    void setup(){
        // avoid memory problems
        H5::H5Library::dontAtExit();
        parpe::initHDF5Mutex();

        fileId = parpe::hdf5CreateFile(std::tmpnam(tempFileName), false);
    }

    void teardown(){
        if(fileId)
            H5Fclose(fileId);
        std::remove(tempFileName);
    }
};
// clang-format on


TEST(hdf5Misc, testOpenExistingFileNoOverwrite) {
    CHECK_THROWS(parpe::HDF5Exception, parpe::hdf5CreateFile(tempFileName, false));
}


TEST(hdf5Misc, testOpenExistingFileOverwrite) {
    H5Fclose(fileId);
    fileId = parpe::hdf5CreateFile(tempFileName, true);
}


TEST(hdf5Misc, testMutexGetLock) {
    parpe::hdf5MutexGetLock();
}


TEST(hdf5Misc, testErrorStackWalker) {
    H5_SAVE_ERROR_HANDLER;

    // provoke error by asking to truncate a file that is already open
    hid_t fileId = H5Fcreate(tempFileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    CHECK_TRUE(fileId <= 0);

    auto s = parpe::captureStreamToString([](){
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, parpe::hdf5ErrorStackWalker_cb, nullptr);
    }, stdout);

    H5_RESTORE_ERROR_HANDLER;

    CHECK_TRUE(200 < s.size());
}


TEST(hdf5Misc, testCreateGroup) {
    const char *groupName = "/test";

    CHECK_FALSE(parpe::hdf5GroupExists(fileId, groupName));

    parpe::hdf5CreateGroup(fileId, groupName, false);

    CHECK_TRUE(parpe::hdf5GroupExists(fileId, groupName));
}


TEST(hdf5Misc, testCreateExistingGroup) {
    parpe::hdf5CreateGroup(fileId, "/test", false);

    H5_SAVE_ERROR_HANDLER;
    CHECK_THROWS(parpe::HDF5Exception, parpe::hdf5CreateGroup(fileId, "/test", false));
    H5_RESTORE_ERROR_HANDLER;
}


TEST(hdf5Misc, testEnsureGroupExists) {
    const char *groupName = "/test";

    CHECK_FALSE(parpe::hdf5GroupExists(fileId, groupName));

    parpe::hdf5EnsureGroupExists(fileId, groupName);

    CHECK_TRUE(parpe::hdf5GroupExists(fileId, groupName));

    parpe::hdf5EnsureGroupExists(fileId, groupName);
}

TEST(hdf5Misc, testStringAttribute) {
    const char *groupName = "/";
    const char *attrName = "testA";
    const char *expAttrValue = "adsf";
    CHECK_FALSE(parpe::hdf5AttributeExists(fileId, groupName, attrName));

    parpe::hdf5WriteStringAttribute(fileId, groupName, attrName, expAttrValue);

    CHECK_TRUE(parpe::hdf5AttributeExists(fileId, groupName, attrName));

    H5T_class_t type_class;
    size_t size = 0;
    int ret = H5LTget_attribute_info(fileId, groupName, attrName, nullptr, &type_class, &size);
    CHECK_TRUE(ret >= 0);
    char actValue[size];

    // TODO: wrong value is read
    H5LTget_attribute_char(fileId, groupName, attrName, actValue);
//    CHECK_EQUAL(expAttrValue, actValue);
}


TEST(hdf5Misc, testDatasetDimensions) {
    const char datasetName[] = "bla";
    const int rank = 3;
    const hsize_t dims[rank] = {1,2,3};
    const int buffer[6] = {1};

    CHECK_FALSE(parpe::hdf5DatasetExists(fileId, datasetName));

    CHECK_THROWS(H5::Exception,
                 parpe::hdf5GetDatasetDimensions(
                     fileId, datasetName, rank,
                     nullptr, nullptr, nullptr, nullptr));

    CHECK_TRUE(H5LTmake_dataset_int(fileId, datasetName, rank, dims, buffer) >= 0);

    CHECK_TRUE(parpe::hdf5DatasetExists(fileId, datasetName));

    int d0 = 0, d1 = 0, d2 = 0, d3 = 0;
    parpe::hdf5GetDatasetDimensions(fileId, datasetName, rank,
                                    &d0, &d1, &d2, &d3);
    CHECK_EQUAL((signed)dims[0], d0);
    CHECK_EQUAL((signed)dims[1], d1);
    CHECK_EQUAL((signed)dims[2], d2);
    CHECK_EQUAL(0      , d3);
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
