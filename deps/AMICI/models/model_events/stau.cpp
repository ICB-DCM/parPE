
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_events{

void stau_model_events(double *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *sx, const int ip, const int ie) {
switch (ip) {
  case 0: {
    switch(ie) { 
        case 0: {
  stau[0] = -(sx[1]-sx[2])/(h[2]+x[2]-p[2]*x[1]+p[1]*x[0]*exp(t*(-1.0/1.0E1))-1.0);

        } break;

        case 1: {
  stau[0] = -(sx[0]-sx[2])/(h[2]+x[2]+p[0]*x[0]*(h[3]-1.0)-1.0);

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

        case 4: {

        } break;

        case 5: {

        } break;

    } 

  } break;

  case 1: {
    switch(ie) { 
        case 0: {
  stau[0] = -(sx[1]-sx[2])/(h[2]+x[2]-p[2]*x[1]+p[1]*x[0]*exp(t*(-1.0/1.0E1))-1.0);

        } break;

        case 1: {
  stau[0] = -(sx[0]-sx[2])/(h[2]+x[2]+p[0]*x[0]*(h[3]-1.0)-1.0);

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

        case 4: {

        } break;

        case 5: {

        } break;

    } 

  } break;

  case 2: {
    switch(ie) { 
        case 0: {
  stau[0] = -(sx[1]-sx[2])/(h[2]+x[2]-p[2]*x[1]+p[1]*x[0]*exp(t*(-1.0/1.0E1))-1.0);

        } break;

        case 1: {
  stau[0] = -(sx[0]-sx[2])/(h[2]+x[2]+p[0]*x[0]*(h[3]-1.0)-1.0);

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

        case 4: {

        } break;

        case 5: {

        } break;

    } 

  } break;

  case 3: {
    switch(ie) { 
        case 0: {
  stau[0] = -(sx[1]-sx[2])/(h[2]+x[2]-p[2]*x[1]+p[1]*x[0]*exp(t*(-1.0/1.0E1))-1.0);

        } break;

        case 1: {
  stau[0] = -(sx[0]-sx[2])/(h[2]+x[2]+p[0]*x[0]*(h[3]-1.0)-1.0);

        } break;

        case 2: {

        } break;

        case 3: {
  stau[0] = 1.0;

        } break;

        case 4: {

        } break;

        case 5: {
  stau[0] = 1.0;

        } break;

    } 

  } break;

}
}

} // namespace model_model_events

} // namespace amici

