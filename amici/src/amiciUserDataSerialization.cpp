#include "amiciUserDataSerialization.h"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>

char *serializeAmiciUserData(const UserData *udata, int *size) {
    std::string serialized;
    boost::iostreams::back_insert_device<std::string> inserter(serialized);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string>>
        s(inserter);
    boost::archive::binary_oarchive oar(s);
    oar << *udata;
    s.flush();

    char *charBuffer = new char[serialized.size()];
    memcpy(charBuffer, serialized.data(), serialized.size());

    if (size)
        *size = serialized.size();

    return charBuffer;
}

UserData deserializeAmiciUserData(const char *buffer, int size) {
    boost::iostreams::basic_array_source<char> device(buffer, size);
    boost::iostreams::stream<boost::iostreams::basic_array_source<char>> s(
        device);
    boost::archive::binary_iarchive iar(s);
    UserData udata;
    iar >> udata;

    return udata;
}
