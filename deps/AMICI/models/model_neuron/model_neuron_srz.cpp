
#include "amici/symbolic_functions.h"
#include "amici/amici_defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void srz_model_neuron(double *srz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip) {
switch (ip) {
  case 0: {
    switch(ie) { 
        case 0: {
  srz[0] = sx[0];

        } break;

    } 

  } break;

  case 1: {
    switch(ie) { 
        case 0: {
  srz[0] = sx[0];

        } break;

    } 

  } break;

  case 2: {
    switch(ie) { 
        case 0: {
  srz[0] = sx[0];

        } break;

    } 

  } break;

  case 3: {
    switch(ie) { 
        case 0: {
  srz[0] = sx[0];

        } break;

    } 

  } break;

}
}

