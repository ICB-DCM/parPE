#ifndef AMICIUSERDATASERIALIZATION_H
#define AMICIUSERDATASERIALIZATION_H

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <udata.h>

char *serializeAmiciUserData(const UserData *udata, int *size);

UserData deserializeAmiciUserData(const char *buffer, int size);

namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive & ar, UserData & u, const unsigned int version)
{
    int size = 0;
    ar & const_cast<int &>(u.np);
    ar & const_cast<int &>(u.nk);
    ar & const_cast<int &>(u.ny);
    ar & const_cast<int &>(u.nytrue);
    ar & const_cast<int &>(u.nx);
    ar & const_cast<int &>(u.nxtrue);
    ar & const_cast<int &>(u.nz);
    ar & const_cast<int &>(u.nztrue);
    ar & const_cast<int &>(u.ne);
    ar & const_cast<int &>(u.nw);
    ar & const_cast<int &>(u.ndwdx);
    ar & const_cast<int &>(u.ndwdp);
    ar & const_cast<int &>(u.nnz);
    ar & const_cast<AMI_o2mode &>(u.o2mode);
    ar & const_cast<int &>(u.ng);
    ar & const_cast<int &>(u.ubw);
    ar & const_cast<int &>(u.lbw);
    ar & u.pscale;
    ar & u.nmaxevent;
    ar & u.nplist;
    ar & u.nt;
    ar & u.tstart;
    ar & u.sensi;
    ar & u.atol;
    ar & u.rtol;
    ar & u.maxsteps;
    ar & u.ism;
    ar & u.sensi_meth;
    ar & u.linsol;
    ar & u.interpType;
    ar & u.lmm;
    ar & u.iter;
    ar & u.stldet;
    ar & u.ordering;
    ar & u.nan_dxdotdp;
    ar & u.nan_J;
    ar & u.nan_JSparse;
    ar & u.nan_xdot;
    ar & u.nan_xBdot;
    ar & u.nan_qBdot;

    if (Archive::is_loading::value) {
        assert(u.qpositivex == NULL);
        ar & size;
        if(size)
            u.qpositivex = new double[size];
    } else {
        size = u.qpositivex == NULL ? 0 : u.nx;
        ar & size;
    }
    ar & make_array<double>(u.qpositivex, size);

    if (Archive::is_loading::value) {
        assert(u.plist == NULL);
        ar & size;
        if(size)
            u.plist = new int[size];
    } else {
        size = u.plist == NULL ? 0 : u.nplist;
        ar & size;
    }
    ar & make_array<int>(u.plist, size);

    if (Archive::is_loading::value) {
        assert(u.p == NULL);
        ar & size;
        if(size)
            u.p = new double[size];
    } else {
        size = u.p == NULL ? 0 : u.np;
        ar & size;
    }
    ar & make_array<double>(u.k, size);

    if (Archive::is_loading::value) {
        assert(u.k == NULL);
        ar & size;
        if(size)
            u.k = new double[size];
    } else {
        size = u.k == NULL ? 0 : u.nk;
        ar & size;
    }

    ar & make_array<double>(u.ts, size);

    if (Archive::is_loading::value) {
        assert(u.ts == NULL);
        ar & size;
        if(size)
            u.ts = new double[size];
    } else {
        size = u.ts == NULL ? 0 : u.nt;
        ar & size;
    }
    ar & make_array<double>(u.ts, size);

    if (Archive::is_loading::value) {
        assert(u.pbar == NULL);
        ar & size;
        if(size)
            u.pbar = new double[size];
    } else {
        size = u.pbar == NULL ? 0 : u.np;
        ar & size;
    }
    ar & make_array<double>(u.pbar, size);

    if (Archive::is_loading::value) {
        assert(u.xbar == NULL);
        ar & size;
        if(size)
            u.xbar = new double[size];
    } else {
        size = u.ts == NULL ? 0 : u.nx;
        ar & size;
    }
    ar & make_array<double>(u.xbar, size);

    if (Archive::is_loading::value) {
        assert(u.x0data == NULL);
        ar & size;
        if(size)
            u.x0data = new double[size];
    } else {
        size = u.x0data == NULL ? 0 : u.nx;
        ar & size;
    }
    ar & make_array<double>(u.x0data, size);

    if (Archive::is_loading::value) {
        assert(u.idlist == NULL);
        ar & size;
        if(size)
            u.idlist = new double[size];
    } else {
        size = u.idlist == NULL ? 0 : u.nx;
        ar & size;
    }
    ar & make_array<double>(u.idlist, size);

    if (Archive::is_loading::value) {
        assert(u.sx0data == NULL);
        ar & size;
        if(size)
            u.sx0data = new double[size];
    } else {
        size = u.sx0data == NULL ? 0 : u.np * u.nx;
        ar & size;
    }
    ar & make_array<double>(u.sx0data, size);

    if (Archive::is_loading::value) {
        assert(u.z2event == NULL);
        ar & size;
        if(size)
            u.z2event = new double[size];
    } else {
        size = u.ts == NULL ? 0 : u.nz;
        ar & size;
    }
    ar & make_array<double>(u.z2event, size);

    if (Archive::is_loading::value) {
        assert(u.h == NULL);
        ar & size;
        if(size)
            u.h = new double[size];
    } else {
        size = u.h == NULL ? 0 : u.ne;
        ar & size;
    }
    ar & make_array<double>(u.h, size);
}

} // namespace serialization
} // namespace boost


#endif // AMICIUSERDATASERIALIZATION_H
