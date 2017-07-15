#include <iostream>
#include <fstream>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

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
    CHECK_EQUAL(u.nxtrue, v.nxtrue);
    CHECK_EQUAL(u.nk, v.nk);
    CHECK_EQUAL(u.ny, v.ny);
    CHECK_EQUAL(u.nytrue, v.nytrue);
    CHECK_EQUAL(u.nx, v.nx);
    CHECK_EQUAL(u.nztrue, v.nztrue);
    CHECK_EQUAL(u.ne, v.ne);
    CHECK_EQUAL(u.nJ, v.nJ);
    CHECK_EQUAL(u.nw, v.nw);
    CHECK_EQUAL(u.ndwdx, v.ndwdx);
    CHECK_EQUAL(u.ndwdp, v.ndwdp);
    CHECK_EQUAL(u.nnz, v.nnz);
    CHECK_EQUAL(u.ubw, v.ubw);
    CHECK_EQUAL(u.lbw, v.lbw);
    CHECK_EQUAL(u.o2mode, v.o2mode);
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
    CHECK_EQUAL(u.nan_dxdotdp, v.nan_dxdotdp);
    CHECK_EQUAL(u.nan_J, v.nan_J);
    CHECK_EQUAL(u.nan_JSparse, v.nan_JSparse);
    CHECK_EQUAL(u.nan_xdot, v.nan_xdot);
    CHECK_EQUAL(u.nan_xBdot, v.nan_xBdot);
    CHECK_EQUAL(u.nan_qBdot, v.nan_qBdot);

}

TEST_GROUP(serialization)
{
    void setup() {
    }

    void teardown() {
    }
};

TEST(serialization, test) {

    UserData u(1,2,3,4,5,6,7,0,0,0,0,0,0,0,0,0,AMICI_SCALING_LN, AMICI_O2MODE_FULL);
//    printUserData(&u);
    {
        std::ofstream ofs( "sstore.dat" );
        boost::archive::text_oarchive oar(ofs);
        oar & u;
    }
    {
        std::ifstream ifs( "sstore.dat" );
        boost::archive::text_iarchive iar(ifs);
        UserData v;
        iar & v;
        checkUserDataEqual(u, v);
//        printUserData(&v);

    }

}

TEST(serialization, test2) {

    UserData u(1,2,3,4,5,6,7,0,0,0,0,0,0,0,0,0,AMICI_SCALING_LN, AMICI_O2MODE_FULL);
//    printUserData(&u);printf("\n");

    std::string serialized;

    {
        boost::iostreams::back_insert_device<std::string> inserter(serialized);
        boost::iostreams::stream<boost::iostreams::back_insert_device<std::string> > s(inserter);
        boost::archive::binary_oarchive oar(s);
        oar << u;
        s.flush();
    }
    {
        boost::iostreams::basic_array_source<char> device(serialized.data(), serialized.size());
        boost::iostreams::stream<boost::iostreams::basic_array_source<char> > s(device);
        boost::archive::binary_iarchive iar(s);
        UserData v;
        iar >> v;
//        printUserData(&v);
        checkUserDataEqual(u, v);
    }
}

TEST(serialization, test3) {

    UserData u(1,2,3,4,5,6,7,0,0,0,0,0,0,0,0,0,AMICI_SCALING_LN, AMICI_O2MODE_FULL);
    u.print();printf("\n");

    int length;
    char *buf = serializeAmiciUserData(&u, &length);
    
    UserData v = deserializeAmiciUserData(buf, length);

    free(buf);
    v.print();
    checkUserDataEqual(u, v);
}



