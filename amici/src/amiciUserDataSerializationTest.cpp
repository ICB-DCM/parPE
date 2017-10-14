#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <fstream>
#include <iostream>

#include <amici_model.h>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include "amiciUserDataSerialization.h"
#include <iomanip>
#include <udata.h>
#include <testingMisc.h>

TEST_GROUP(userDataSerialization){void setup(){}

                          void teardown(){}};

TEST(userDataSerialization, testFile) {

    UserData u(1, 2, 3);
    {
        std::ofstream ofs("sstore.dat");
        boost::archive::text_oarchive oar(ofs);
        oar &u;
    }
    {
        std::ifstream ifs("sstore.dat");
        boost::archive::text_iarchive iar(ifs);
        UserData v;
        iar &v;
        checkUserDataEqual(u, v);
    }
}

TEST(userDataSerialization, testString) {

    UserData u(1, 2, 3);

    std::string serialized = serializeToString(u);

    checkUserDataEqual(u, deserializeFromString<UserData>(serialized));
}

TEST(userDataSerialization, testChar) {

    UserData u(1, 2, 3);
    u.p = new double[2];
    u.p[0] = 1;
    u.p[1] = 2;

    int length;
    char *buf = serializeToChar(&u, &length);

    UserData v = deserializeFromChar<UserData>(buf, length);

    delete[] buf;
    checkUserDataEqual(u, v);
}

TEST_GROUP(returnDataSerialization){void setup(){}

                          void teardown(){}};

TEST(returnDataSerialization, testString) {
    UserData u(1, 2, 3);
    Model m(1,2,3,3,4,5,6,7,8,9,10, 11, 12, 13,14,15, AMICI_O2MODE_NONE);


    ReturnData r(&u, &m);

    std::string serialized = serializeToString(u);

    checkReturnDataEqual(r, deserializeFromString<ReturnData>(serialized));
}

