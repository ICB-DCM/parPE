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

TEST(serialization, test2) {

    UserData u(1,2,3,4,0,0,0,0,0,0,0,0,0,0,0,0,AMI_SCALING_LN, AMI_O2MODE_FULL);
    printUserData(&u);

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
        printUserData(&v);
    }
}


