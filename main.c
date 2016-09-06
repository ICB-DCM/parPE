#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#include "wrapfunctions.h"
#include "include/symbolic_functions.h"
#include "include/amici.h"
#include "src/ami_hdf5.h"

#define EXPERIMENT_INDEX_CONTROL -1
#define XDOT_REL_TOLERANCE 1e-6
#define NUM_FIXED_PARAMS 765
#define NUM_STATE_VARIABLES 1230
#define NUM_CELL_LINES 56
#define NUM_EXPERIMENTS_PER_CELLLINE 9

UserData getMyUserData() {
    // common udata for all?

    const char* fileName = "/home/dweindl/src/CanPathProSSH/dw/test.h5";;
    UserData udata; /* User udata structure */
    udata = (UserData) malloc(sizeof *udata);
    if (udata == NULL)
        return(NULL);

    init_modeldims(udata);

    hid_t file_id = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
    //status = H5Aread (attr_id, mem_type_id, buf);

    const char* optionsObject = "/options";
    atol = getDoubleScalarAttribute(file_id, optionsObject, "atol");
    rtol = getDoubleScalarAttribute(file_id, optionsObject, "rtol");
    maxsteps = getIntScalarAttribute(file_id, optionsObject, "maxsteps");
    tstart = getDoubleScalarAttribute(file_id, optionsObject, "tstart");
    lmm = getIntScalarAttribute(file_id, optionsObject, "lmm");
    iter = getIntScalarAttribute(file_id, optionsObject, "iter");
    linsol = getIntScalarAttribute(file_id, optionsObject, "linsol");
    stldet = getIntScalarAttribute(file_id, optionsObject, "stldet");
    interpType = getIntScalarAttribute(file_id, optionsObject, "interpType");
    //lmmB = getDoubleScalarAttribute(file_id, optionsObject, "lmmB");
    //iterB = getDoubleScalarAttribute(file_id, optionsObject, "iterB");
    ism = getIntScalarAttribute(file_id, optionsObject, "ism");
    sensi_meth = getIntScalarAttribute(file_id, optionsObject, "sensi_meth");
    sensi = getIntScalarAttribute(file_id, optionsObject, "sensi");
    nmaxevent = getIntScalarAttribute(file_id, optionsObject, "nmaxevent");
    ordering = getIntScalarAttribute(file_id, optionsObject, "ordering");
    //ss = getDoubleScalarAttribute(file_id, optionsObject, "ss");

    hsize_t length;
    getDoubleArrayAttribute(file_id, optionsObject, "qpositivex", &qpositivex, &length);
    // getIntArrayAttribute(file_id, optionsObject, "sens_ind", sens_ind, &tmp);
    //    'sx0':  []

    const char* dataObject = "/data";
    getDoubleArrayAttribute(file_id, dataObject, "theta", &p, &length);
    np = length;

    getDoubleArrayAttribute(file_id, dataObject, "kappa", &k, &length);
    assert(length == getIntScalarAttribute(file_id, dataObject, "nk"));
    getDoubleArrayAttribute(file_id, dataObject, "ts", &ts, &length);
    nt = getIntScalarAttribute(file_id, dataObject, "nt");
    assert(length == nt);
    assert(ny == getIntScalarAttribute(file_id, dataObject, "ny"));
    assert(nz == getIntScalarAttribute(file_id, dataObject, "nz"));

    k = malloc(sizeof(double) * (NUM_FIXED_PARAMS + NUM_STATE_VARIABLES));

    /* plist, matlab: fifth argument */
    // parameter ordering
    plist = malloc(np * sizeof(int));
    for (int i = 0; i < np; i++) {
        plist[i] = i;
    }

    /* Options ; matlab: fourth argument   */
    z2event = malloc(sizeof(realtype) * ne);
    for(int i = 0; i < ne; ++i)
        z2event[i] = i;

    idlist = malloc(sizeof(realtype) * np);
    for(int i = 0; i < np; ++i)
        idlist[i] = 0;
    //user-provided sensitivity initialisation. this should be a matrix of dimension [#states x #parameters] default is sensitivity initialisation based on the derivative of the state initialisation
    b_sx0 = FALSE;
    sx0data = 0;

    /* pbarm parameterscales ; matlab: sixth argument*/
    pbar = malloc(sizeof(realtype) * np);
    ones(pbar, np);

    //    /* xscale, matlab: seventh argument */
    //    xbar = mxGetPr(prhs[6]);
    xbar = 0; // xscale

//    double *am_qpositivex;
//    /** parameter reordering */
//    int    *am_plist;
//    /** number of parameters */
//    int    am_np;
//    /** number of observables */
//    int    am_ny;
//    /** number of observables in the unaugmented system */
//    int    am_nytrue;
//    /** number of states */
//    int    am_nx;
//    /** number of event outputs */
//    int    am_nz;
//    /** number of event outputs in the unaugmented system */
//    int    am_nztrue;
//    /** number of events */
//    int    am_ne;
//    /** number of timepoints */
//    int    am_nt;
//    /** number of common expressions */
//    int    am_nw;
//    /** number of derivatives of common expressions wrt x */
//    int    am_ndwdx;
//    /** number of derivatives of common expressions wrt p */
//    int    am_ndwdp;
//    /** number of nonzero entries in jacobian */
//    int    am_nnz;

//    /** scaling of parameters */
//    double *am_pbar;
//    /** scaling of states */
//    double *am_xbar;
//    /** flag array for DAE equations */
//    double *am_idlist;

//    /** upper bandwith of the jacobian */
//    int am_ubw;
//    /** lower bandwith of the jacobian */
//    int am_lbw;

//    /** flag for sensitivity initialisation */
//    /*!
//     * flag which determines whether analytic sensitivities initialisation or provided initialisation should be used
//     */
//    booleantype am_bsx0;

//    /** sensitivity initialisation */
//    double *am_sx0data;

//    /** index indicating to which event an event output belongs */
//    double *am_z2event;

//    /** flag indicating whether a certain heaviside function should be active or not */
//    double *am_h;

//    /** tempory storage of Jacobian data across functions */
//    SlsMat am_J;
//    /** tempory storage of dxdotdp data across functions */
//    realtype *am_dxdotdp;
//    /** tempory storage of w data across functions */
//    realtype *am_w;
//    /** tempory storage of dwdx data across functions */
//    realtype *am_dwdx;
//    /** tempory storage of dwdp data across functions */
//    realtype *am_dwdp;
//    /** tempory storage of M data across functions */
//    realtype *am_M;
//    /** tempory storage of dfdx data across functions */
//    realtype *am_dfdx;
//    /** tempory storage of stau data across functions */
//    realtype *am_stau;

//    /** flag indicating whether a NaN in dxdotdp has been reported */
//    booleantype am_nan_dxdotdp;
//    /** flag indicating whether a NaN in J has been reported */
//    booleantype am_nan_J;
//    /** flag indicating whether a NaN in JSparse has been reported */
//    booleantype am_nan_JSparse;
//    /** flag indicating whether a NaN in xdot has been reported */
//    booleantype am_nan_xdot;
//    /** flag indicating whether a NaN in xBdot has been reported */
//    booleantype am_nan_xBdot;
//    /** flag indicating whether a NaN in qBdot has been reported */
//    booleantype am_nan_qBdot;

    H5Fclose(file_id);
    processUserData(udata);

    return(udata);
}


void readFixedParameters(int dataMatrixRow, hid_t file_id, UserData udata) {
    // set kappa

    int rank = 2;
    hsize_t dims[] = {1, NUM_FIXED_PARAMS + NUM_STATE_VARIABLES};
    hid_t memspaceId = H5Screate_simple(rank, dims, NULL);
    hid_t datasetId = H5Dopen2(file_id, "/data/C", H5P_DEFAULT);
    hid_t dataspaceId = H5Dget_space(datasetId);
    hsize_t offset[] = {dataMatrixRow, 0};
    hsize_t count[] = {1, NUM_FIXED_PARAMS + NUM_STATE_VARIABLES};

    //buffer = malloc(sizeof(double) * dims[0] * dims[1]);
    herr_t status = H5Sselect_hyperslab(dataspaceId, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dread(datasetId, H5T_NATIVE_DOUBLE, memspaceId, dataspaceId, H5P_DEFAULT, k);

}

ExpData getExperimentalDataForExperiment(int cellLineIdx, int expIdx, UserData udata) {
    const char *hdffile = "/home/dweindl/src/CanPathProSSH/dw/data.h5";

    ExpData edata = (ExpData) malloc(sizeof *edata);
    if (edata == NULL) {
        return(NULL);
    }

    hid_t file_id = H5Fopen(hdffile, H5F_ACC_RDONLY, H5P_DEFAULT);

    int dataMatrixIdx = -1;

    char path[50];
    int _k = 1;
    assert(_k < 10 && _k >= 0 && cellLineIdx >= 0 && cellLineIdx < 100);
    sprintf(path, "/data/I%d/%d", _k , cellLineIdx);

    hsize_t len;
    double *buffer;

    my = malloc(sizeof(double));
    ysigma = malloc(sizeof(double));
    mz = 0;
    zsigma = 0;

    if(expIdx == EXPERIMENT_INDEX_CONTROL) {
        *my = 1;
        *ysigma = NAN; // or 1?
        int *intBuf;
        getIntArrayAttribute(file_id, path, "ref", &intBuf, &len);
        dataMatrixIdx = intBuf[0];
        free(intBuf);
    } else {
        getDoubleArrayAttribute(file_id, path, "D.Y", &buffer, &len);
        assert(expIdx >= 0 && expIdx < len);
        *my = buffer[expIdx];
        free(buffer);
        getDoubleArrayAttribute(file_id, path, "D.Sigma_Y", &buffer, &len);
        assert(expIdx >= 0 && expIdx < len);
        *ysigma = buffer[expIdx]; // NAN?
        free(buffer);

        int *intBuf;
        getIntArrayAttribute(file_id, path, "exp", &intBuf, &len);
        assert(expIdx >= 0 && expIdx < len);
        dataMatrixIdx = intBuf[0];
        free(intBuf);
    };

    readFixedParameters(dataMatrixIdx, file_id, udata);

    H5Fclose(file_id);

    return(edata);
}

void error(const char *message) { // exit?
    printf("ERROR: %s\n", message);
}

void warning(const char *message) {
    printf("WARNING: %s\n", message);
}

bool reachedSteadystate(double *_xdot, double *_x, int numTimepoints, int numStates, double tolerance) {
    for(int state = 0; state < numStates; ++state) {
        for(int time = 0; time < numTimepoints; ++time) {
            double sensitivity = _xdot[state] / (_x[state] + tolerance); //TODO check dimensions
            if(sensitivity > tolerance) {
                return FALSE;
            }
        }
    }
    return TRUE;
}

void updateInitialConditions(double *oldInitialConditions, const double *newInitialConditions, int count) {
    memcpy(newInitialConditions, oldInitialConditions, count * sizeof(double));
}

ReturnData getSteadyStateSolution(UserData udata, ExpData edata, int *status) {
    ReturnData rdata;
    while (TRUE) {
        rdata = getSimulationResults(udata, edata, status);

        if(*status < 0) {
            error("Failed to integrate."); // TODO add dataset info, case/control, celline
        }

        bool hasSteadystate = reachedSteadystate(xdotdata, xdata, nt, nx, XDOT_REL_TOLERANCE);

        if(hasSteadystate)
            break;

        // use previous solution as initial conditions
        updateInitialConditions(xdata, &k[NUM_FIXED_PARAMS], NUM_STATE_VARIABLES);
    }

    return rdata;
}

double evaluateObjectiveFunction(double* timepoints, double *theta, double lenTheta, int numCellline, int numExperiments) { // (tin,theta,C,I)
    double *logLikelihoodJ = malloc(sizeof(double) * numCellline);
    zeros(logLikelihoodJ, numCellline);

    double *dLogLikelihoodJ = malloc(sizeof(double) * numCellline * lenTheta);
    zeros(dLogLikelihoodJ, numCellline * lenTheta);

    // set options ... as argument
    // init amidata ...

    UserData udata = getMyUserData();
    int status = 0;

    for(int celllineIdx = 1; celllineIdx <= numCellline; ++celllineIdx) {
        printf("Cell line idx %d\n", celllineIdx);
        fflush(stdout);
        ReturnData rdata;
        ExpData edata = getExperimentalDataForExperiment(celllineIdx, EXPERIMENT_INDEX_CONTROL, udata);
        // print & check data


        ReturnData rdataControl = getSteadyStateSolution(udata, edata, &status);
        rdata = rdataControl;
        freeExpData(edata); // reuse?

        // 'case' experiments
        for(int experimentIdx = 0; experimentIdx < numExperiments; ++experimentIdx) {
            printf("\tExperiment idx %d\n", experimentIdx);
            fflush(stdout);

            edata = getExperimentalDataForExperiment(celllineIdx, experimentIdx, udata);
            // use control solution as starting values
            updateInitialConditions(&k[NUM_FIXED_PARAMS], xdata, NUM_STATE_VARIABLES);
            ReturnData rdataCase = getSteadyStateSolution(udata, edata, &status);

            double growthInhibition = *rdataCase->am_llhdata / *rdataControl->am_llhdata;

            // TODO store

            double weightedError = ((growthInhibition - my[0]) / ysigma[0]);
            logLikelihoodJ[celllineIdx] -= 0.5 * weightedError * weightedError;

            free(edata);

//            if(nderiv>1) {
//                dinhib = (sol_drug.sllh.*sol_ref.llh - sol_drug.llh.*sol_ref.sllh)./(sol_ref.llh.^2);
//                dlogLj(:,icl) = dlogLj(:,icl) - ( ( inhib - Y )*( dinhib ) / Sigma_Y^2 );
//            }

        }
    }

    // free udata, edata, rdata

    if(status == 0) {
        double logLikelihood = sum(logLikelihoodJ, numCellline);
//        if(nderiv > 1) {
//            //                 dlogL = sum(dlogLj,2);
//        }
        return logLikelihood;
    }

    return INFINITY;
}

void optimizeModel() {
    bool converged = 1;
    do {
        double timepoints [] = {};
        double theta [] = {};
        double j = evaluateObjectiveFunction(timepoints, theta,
                                             NUM_FIXED_PARAMS + NUM_STATE_VARIABLES,
                                             NUM_CELL_LINES, NUM_EXPERIMENTS_PER_CELLLINE);
        printf("Objective function value %e\n", j);
    } while(!converged);
}

int main(int argc, char **argv)
{
    optimizeModel();
    printf("Test\n");
}


