#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <fstream>
#include <iostream>

#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include "amiciUserDataSerialization.h"
#include <iomanip>
#include <udata.h>

void checkUserDataEqual(UserData u, UserData v) {
    CHECK_EQUAL(u.np, v.np);
    CHECK_EQUAL(u.nx, v.nx);
    CHECK_EQUAL(u.nk, v.nk);
    CHECK_EQUAL(u.nx, v.nx);
    CHECK_EQUAL(u.pscale, v.pscale);

    CHECK_EQUAL(u.nmaxevent, v.nmaxevent);
    CHECK_EQUAL(u.nplist, v.nplist);
    CHECK_EQUAL(u.nt, v.nt);
    CHECK_EQUAL(u.tstart, v.tstart);
    CHECK_EQUAL(u.sensi, v.sensi);
    CHECK_EQUAL(u.atol, v.atol);
    CHECK_EQUAL(u.rtol, v.rtol);
    CHECK_EQUAL(u.maxsteps, v.maxsteps);
    CHECK_EQUAL(u.ism, v.ism);
    CHECK_EQUAL(u.sensi_meth, v.sensi_meth);
    CHECK_EQUAL(u.linsol, v.linsol);
    CHECK_EQUAL(u.interpType, v.interpType);
    CHECK_EQUAL(u.lmm, v.lmm);
    CHECK_EQUAL(u.iter, v.iter);
    CHECK_EQUAL(u.stldet, v.stldet);
    CHECK_EQUAL(u.ordering, v.ordering);
}

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

    std::string serialized;

    {
        boost::iostreams::back_insert_device<std::string> inserter(serialized);
        boost::iostreams::stream<
            boost::iostreams::back_insert_device<std::string>>
            s(inserter);
        boost::archive::binary_oarchive oar(s);
        oar << u;
        s.flush();
    }
    {
        boost::iostreams::basic_array_source<char> device(serialized.data(),
                                                          serialized.size());
        boost::iostreams::stream<boost::iostreams::basic_array_source<char>> s(
            device);
        boost::archive::binary_iarchive iar(s);
        UserData v;
        iar >> v;
        checkUserDataEqual(u, v);
    }
}

TEST(userDataSerialization, testChar) {

    UserData u(1, 2, 3);

    int length;
    char *buf = serializeToChar(&u, &length);

    UserData v = deserializeFromChar<UserData>(buf, length);

    delete[] buf;
    checkUserDataEqual(u, v);
}
