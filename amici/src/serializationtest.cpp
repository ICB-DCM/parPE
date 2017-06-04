#include <iostream>
#include <fstream>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include "amiciUserDataSerialization.h"
#include <iomanip>
#include <udata.h>


TEST_GROUP(serialization)
{
    void setup() {
    }

    void teardown() {
    }
};

TEST(serialization, test) {

    UserData u(1,2,3,4,0,0,0,0,0,0,0,0,0,0,0,0,AMI_SCALING_LN, AMI_O2MODE_FULL);
    printUserData(&u);
    {
        std::ofstream ofs( "store.dat" );
        boost::archive::text_oarchive oar(ofs);
        oar & u;
    }
    {
        std::ifstream ifs( "store.dat" );
        boost::archive::text_iarchive iar(ifs);
        UserData v;
        iar & v;
        printUserData(&v);
    }
}


