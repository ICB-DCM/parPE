#include <math.h>
#include <limits.h>
#include <assert.h>
#include <pthread.h>
#include <execinfo.h>
#include <unistd.h>

#include "dataprovider.h"

#include <wrapfunctions.h>
#include <src/ami_hdf5.h>
#include <include/udata_accessors.h>
#include <include/edata_accessors.h>


#define H5_SAVE_ERROR_HANDLER   herr_t (*old_func)(void*); \
                                void *old_client_data; \
                                H5Eget_auto1(&old_func, &old_client_data); \
                                H5Eset_auto1(NULL, NULL);

#define H5_RESTORE_ERROR_HANDLER H5Eset_auto1(old_func, old_client_data);

static hid_t file_id;

// TODO check docs if this can be avoided
extern pthread_mutex_t mutexHDF;


static void readFixedParameters(int dataMatrixRow, UserData *udata);

// malloc version of ami_hdf5.cpp
static void getDoubleArrayAttributeC(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *length);
static void getIntArrayAttributeC(hid_t file_id, const char* optionsObject, const char* attributeName, int **destination, hsize_t *length);



int initDataProvider(const char *hdf5Filename) {

#if DATA_PROVIDER_H_VERBOSE >= 2
    printf("Entering dataprovider.c initDataProvider '%s'.\n", hdf5Filename); fflush(stdout);
#endif

    H5_SAVE_ERROR_HANDLER;
    file_id = H5Fopen(hdf5Filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    if(file_id < 0) {
        logmessage(LOGLVL_CRITICAL, "initDataProvider failed to open HDF5 file '%s'.", hdf5Filename);
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, NULL);
    }
    H5_RESTORE_ERROR_HANDLER;

    return file_id < 0;
}

void closeDataProvider() {
    H5_SAVE_ERROR_HANDLER
    herr_t status = H5Fclose(file_id);

    if(status< 0) {
        error("closeDataProvider failed to close HDF5 file.");
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, NULL);
    }
    H5_RESTORE_ERROR_HANDLER
}

int getNumMultiStartRuns() {
    // TODO
    return 1;
}

int getNumLocalOptimizationsForMultiStartRun(int multiStartRun) {
    // TODO
    return 2;
}

int getMaxIter() {
    // TODO
    return 1;
}

int getNumGenotypes(datapath path)
{
    // TODO
    return 1;
    return 96;
}

int getExperimentCountForCellline(datapath dpath)
{
    return 1; // TODO testing

    char path[50];
    int _k = 1;
    assert(_k < 10 && _k >= 0 && dpath.idxGenotype >= 0 && dpath.idxGenotype < 100);
    sprintf(path, "/data/I%d/%d", _k , dpath.idxGenotype);

    hsize_t len;
    H5T_class_t type_class;
    size_t type_size;
    herr_t status;

    pthread_mutex_lock(&mutexHDF);
    status = H5LTget_attribute_info(file_id, path, "D.Y", &len, &type_class, &type_size);
    pthread_mutex_unlock(&mutexHDF);

    if(status < 0) {
        fprintf(stderr, "Error in getExperimentCountForCellline: Cannot read attribute '%s' of '%s'\n", "D.Y", path);
        void *array[10];
        size_t size;
        size = backtrace(array, 10);
        backtrace_symbols_fd(array, size, STDERR_FILENO);
    }

    return(len);
}

int getNumberOfSimulationForObjectiveFunction(datapath path) {
    int count = 0;

    int numGenotypes = getNumGenotypes(path);

    for(int genotypeIdx = 0; genotypeIdx < numGenotypes; ++genotypeIdx) {
        path.idxGenotype = genotypeIdx + 1; // starting from 1

        int numExperiments = getExperimentCountForCellline(path);
        count += numExperiments;
    }

    return count;
}

UserData *getMyUserData() {

#if DATA_PROVIDER_H_VERBOSE >= 2
    printf("Entering dataprovider.c getMyUserData.\n"); fflush(stdout);
#endif

    // common udata for all?
    // No, different protein configuration (per cellline) and exppression data (per condition)
    // Rest should be same

    UserData *udata = malloc(sizeof(UserData));
    if (udata == NULL)
        return(NULL);

    init_modeldims(udata);

    pthread_mutex_lock(&mutexHDF);

    H5_SAVE_ERROR_HANDLER

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
    ism = getIntScalarAttribute(file_id, optionsObject, "ism");
    sensi_meth = getIntScalarAttribute(file_id, optionsObject, "sensi_meth");
    sensi = getIntScalarAttribute(file_id, optionsObject, "sensi");
    nmaxevent = getIntScalarAttribute(file_id, optionsObject, "nmaxevent");
    ordering = getIntScalarAttribute(file_id, optionsObject, "ordering");

    hsize_t length;
    getDoubleArrayAttributeC(file_id, optionsObject, "qpositivex", &qpositivex, &length);

    const char* dataObject = "/modelinfo";

    getDoubleArrayAttributeC(file_id, dataObject, "theta", &p, &length);
    np = length;

    getDoubleArrayAttributeC(file_id, dataObject, "kappa", &k, &length);
    // assert(length == getIntScalarAttribute(file_id, dataObject, "nk")); // nk is only written during first simulation run
    getDoubleArrayAttributeC(file_id, dataObject, "ts", &ts, &length);
    nt = getIntScalarAttribute(file_id, dataObject, "nt");
    assert(length == nt);

    assert(ny == getIntScalarAttribute(file_id, dataObject, "ny"));
    assert(nz == getIntScalarAttribute(file_id, dataObject, "nz"));


    if(H5Eget_num(H5E_DEFAULT)) {
        error("Problem with reading udata\n");
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, NULL);
    }
    H5_RESTORE_ERROR_HANDLER

    pthread_mutex_unlock(&mutexHDF);

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

    // TODO b_x0  x0data
    b_x0 = FALSE; // just use x0data as pointer to k, set b_x0 = 0 to not try to deallocate x0data memory
    x0data = &k[NUM_FIXED_PARAMS];

    b_sx0 = FALSE;
    sx0data = 0;

    pbar = malloc(sizeof(realtype) * np);
    ones(pbar, np);

    xbar = 0; // xscale
    //    'sx0':  []  see AMI_setup
    //    int    am_nt;
    h = 0;
    udata->am_dxdotdp = 0;
    udata->am_dwdx = 0;
    udata->am_dwdp = 0;
    udata->am_M = 0;
    udata->am_dfdx = 0;
    udata->am_stau = 0;
    stau_tmp = 0;

    udata->am_pscale = AMI_SCALING_LOG10;

// otherwise not able to free amici-allocated space... TODO need to refactor
//    processUserData(udata);
//    void processUserData(UserData *udata) {
        if (nx>0) {
            /* initialise temporary jacobian storage */
            tmp_J = SparseNewMat(nx,nx,nnz,CSC_MAT);
            M_tmp = malloc(sizeof(realtype) * nx * nx);
            dfdx_tmp = malloc(sizeof(realtype) * nx*nx);
        }
        if (sensi>0) {
            /* initialise temporary dxdotdp storage */
            tmp_dxdotdp = malloc(sizeof(realtype) * nx*np);
        }
        if (ne>0) {
            /* initialise temporary stau storage */
            stau_tmp = malloc(sizeof(realtype) * np);
        }


        w_tmp = malloc(sizeof(realtype) * nw);
        dwdx_tmp = malloc(sizeof(realtype) * ndwdx);
        dwdp_tmp = malloc(sizeof(realtype) * ndwdp);
//    }

    return(udata);
}

void readFixedParameters(int dataMatrixRow, UserData *udata) {

#if DATA_PROVIDER_H_VERBOSE >= 2
    printf("Entering readFixedParameters.\n"); fflush(stdout);
#endif

    // set kappa
    int rank = 2;
    hsize_t dims[] = {1, NUM_FIXED_PARAMS + NUM_STATE_VARIABLES};

    pthread_mutex_lock(&mutexHDF);

    H5_SAVE_ERROR_HANDLER

    hid_t memspaceId = H5Screate_simple(rank, dims, NULL);
    hid_t datasetId = H5Dopen2(file_id, "/data/C", H5P_DEFAULT);
    hid_t dataspaceId = H5Dget_space(datasetId);
    hsize_t offset[] = {dataMatrixRow, 0};
    hsize_t count[] = {1, NUM_FIXED_PARAMS + NUM_STATE_VARIABLES};

    H5Sselect_hyperslab(dataspaceId, H5S_SELECT_SET, offset, NULL, count, NULL);
    H5Dread(datasetId, H5T_NATIVE_DOUBLE, memspaceId, dataspaceId, H5P_DEFAULT, k);

    if(H5Eget_num(H5E_DEFAULT)) {
        error("Problem in readFixedParameters\n");
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, NULL);
        abort();
    }

    H5_RESTORE_ERROR_HANDLER

    pthread_mutex_unlock(&mutexHDF);
}

ExpData *getExperimentalDataForExperiment(datapath dpath, UserData *udata) {

#if DATA_PROVIDER_H_VERBOSE >= 2
    printf("Entering getExperimentalDataForExperiment "); printDatapath(dpath); fflush(stdout);
#endif

    char path[50];
    int _k = 1; // TODO dpath.idx...
    assert(_k < 10 && _k >= 0);
    assert(dpath.idxGenotype > 0 && dpath.idxGenotype < 100); // to fit into 'path'
    sprintf(path, "/data/I%d/%d", _k , dpath.idxGenotype);

    ExpData *edata = malloc(sizeof(ExpData));
    my = malloc(sizeof(double)); // TODO larger Y
    ysigma = malloc(sizeof(double)); // TODO larger Y
    mz = 0;
    zsigma = 0;

    // find index in "C" matrix
    int dataMatrixIdx = -1;

    hsize_t len;
    double *buffer;
    pthread_mutex_lock(&mutexHDF);

    if(dpath.idxExperiment == EXPERIMENT_INDEX_CONTROL) {
        my[0] = 1;
        ysigma[0] = NAN; // TODO or 1?
        int *intBuf;
        getIntArrayAttributeC(file_id, path, "ref", &intBuf, &len);
        dataMatrixIdx = intBuf[0];
        free(intBuf);
    } else {
        getDoubleArrayAttributeC(file_id, path, "D.Y", &buffer, &len);
        assert(dpath.idxExperiment >= 0 && dpath.idxExperiment < len);
        *my = buffer[dpath.idxExperiment];
        free(buffer);
        getDoubleArrayAttributeC(file_id, path, "D.Sigma_Y", &buffer, &len);
        assert(dpath.idxExperiment >= 0 && dpath.idxExperiment < len);
        *ysigma = buffer[dpath.idxExperiment]; // NAN?
        free(buffer);

        int *intBuf;
        getIntArrayAttributeC(file_id, path, "exp", &intBuf, &len);
        assert(dpath.idxExperiment >= 0 && dpath.idxExperiment < len);
        dataMatrixIdx = intBuf[0];
        free(intBuf);
    };

    b_expdata = TRUE;

    readFixedParameters(dataMatrixIdx, udata);
    pthread_mutex_unlock(&mutexHDF);

    return edata;
}

void myFreeExpData(ExpData *edata) {
    if(edata) {
        if(my) free(my);
        if(ysigma) free(ysigma);
        if(mz) free(mz);
        if(zsigma) free(zsigma);
        free(edata);
    }
}


void freeUserDataC(UserData *udata) {
    if(udata) {
        if(qpositivex) free(qpositivex);
        if(p) free(p);
        if(k) free(k);
        if(ts) free(ts);
        if(pbar) free(pbar);
        if(xbar) free(xbar);
        if(idlist) free(idlist);
        if(sx0data) free(sx0data);
        if(z2event) free(z2event);
        if(plist) free(plist);
        if(h) free(h);
        if(tmp_dxdotdp) free(tmp_dxdotdp);
        if(w_tmp) free(w_tmp);
        if(dwdx_tmp) free(dwdx_tmp);
        if(dwdp_tmp) free(dwdp_tmp);
        if(M_tmp) free(M_tmp);
        if(dfdx_tmp) free(dfdx_tmp);
        if(stau_tmp) free(stau_tmp);
        SparseDestroyMat(tmp_J);

        free(udata);
    }
}

void getThetaLowerBounds(datapath dataPath, double *buffer, AMI_parameter_scaling scaling)
{
    if(scaling == AMI_SCALING_LOG10)
        fillArray(buffer, getLenTheta(), -2);
    else
        fillArray(buffer, getLenTheta(), 1e-2);
}

void getThetaUpperBounds(datapath dataPath, double *buffer, AMI_parameter_scaling scaling)
{
    if(scaling == AMI_SCALING_LOG10)
        fillArray(buffer, getLenTheta(), 2);
    else
        fillArray(buffer, getLenTheta(), 1e2);
}

void getFeasibleInitialThetaFromFile(datapath dataPath, double *buffer, AMI_parameter_scaling scaling) {
    char path[100];
    sprintf(path, "/crossvalidations/%d/multistarts/%d/initialTheta", dataPath.idxMultiStart, dataPath.idxLocalOptimization  % 5); // TODO %5 to reuse starting points for scaling tests

    logmessage(LOGLVL_INFO, "Reading initial theta from %s", path);

    pthread_mutex_lock(&mutexHDF);
    H5LTread_dataset_double(file_id, path, buffer);
    pthread_mutex_unlock(&mutexHDF);
}

void getRandomInitialThetaFromFile(datapath dataPath, double *buffer, AMI_parameter_scaling scaling) {
    // TODO matlab writes fortran style, so this reads rowwise, transpose...

    logmessage(LOGLVL_INFO, "Reading random initial theta %d from /randomstarts", dataPath.idxLocalOptimization);

    int rank = 2;
    hsize_t dims[] = {1, NUM_OPTIMIZATION_PARAMS};

    pthread_mutex_lock(&mutexHDF);
    H5_SAVE_ERROR_HANDLER

    hid_t memspaceId = H5Screate_simple(rank, dims, NULL);
    hid_t datasetId = H5Dopen2(file_id, "/randomstarts", H5P_DEFAULT);
    hid_t dataspaceId = H5Dget_space(datasetId);
    hsize_t offset[] = {dataPath.idxLocalOptimization, 0};
    hsize_t count[] = {1, NUM_OPTIMIZATION_PARAMS};

    H5Sselect_hyperslab(dataspaceId, H5S_SELECT_SET, offset, NULL, count, NULL);
    H5Dread(datasetId, H5T_NATIVE_DOUBLE, memspaceId, dataspaceId, H5P_DEFAULT, buffer);

    if(H5Eget_num(H5E_DEFAULT)) {
        error("Problem in getRandomInitialThetaFromFile\n");
        H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, hdf5ErrorStackWalker_cb, NULL);
    }
    abort();

    H5_RESTORE_ERROR_HANDLER
    pthread_mutex_unlock(&mutexHDF);
}

int getLenTheta()
{
    return NUM_OPTIMIZATION_PARAMS;
}

int getLenKappa()
{
    return NUM_FIXED_PARAMS;
}

void printDatapath(datapath path)
{
    printf("%d.%d.%d.%d.%d\n", path.idxMultiStart, path.idxLocalOptimization, path.idxLocalOptimizationIteration, path.idxGenotype, path.idxExperiment);
}

void sprintDatapath(char *buffer, datapath path)
{
    sprintf(buffer, "%d.%d.%d.%d.%d", path.idxMultiStart, path.idxLocalOptimization, path.idxLocalOptimizationIteration, path.idxGenotype, path.idxExperiment);
}

herr_t hdf5ErrorStackWalker_cb(unsigned int n, const H5E_error_t *err_desc, void *client_data)
{
    assert (err_desc);
    const int		indent = 2;

    const char *maj_str = H5Eget_major(err_desc->maj_num);
    const char *min_str = H5Eget_minor(err_desc->min_num);

    logmessage(LOGLVL_CRITICAL, "%*s#%03d: %s line %u in %s(): %s",
         indent, "", n, err_desc->file_name, err_desc->line,
         err_desc->func_name, err_desc->desc);
    logmessage(LOGLVL_CRITICAL, "%*smajor(%02d): %s",
         indent*2, "", err_desc->maj_num, maj_str);
    logmessage(LOGLVL_CRITICAL, "%*sminor(%02d): %s",
         indent*2, "", err_desc->min_num, min_str);

    return 0;
}

void getDoubleArrayAttributeC(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *length) {
    H5T_class_t type_class;
    size_t type_size;
    herr_t status;

    pthread_mutex_lock(&mutexHDF);

    status = H5LTget_attribute_info(file_id, optionsObject, attributeName, length, &type_class, &type_size);
    if(status < 0) {
        fprintf(stderr, "Error in getDoubleArrayAttributeC: Cannot read attribute '%s' of '%s'\n", attributeName, optionsObject);
        void *array[10];
        size_t size;
        size = backtrace(array, 10);
        backtrace_symbols_fd(array, size, STDERR_FILENO);
    }

    *destination = malloc(sizeof(double) * (*length)); // vs. type_size
    status = H5LTget_attribute_double(file_id, optionsObject, attributeName, *destination);
    if(status < 0)
        fprintf(stderr, "Error in getDoubleArrayAttributeC: Cannot read attribute '%s' of '%s'\n", attributeName, optionsObject);

    pthread_mutex_unlock(&mutexHDF);
}

void getIntArrayAttributeC(hid_t file_id, const char* optionsObject, const char* attributeName, int **destination, hsize_t *length) {
    H5T_class_t type_class;
    size_t type_size;

    pthread_mutex_lock(&mutexHDF);

    H5LTget_attribute_info(file_id, optionsObject, attributeName, length, &type_class, &type_size);
#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %d: ", attributeName, *length);
#endif
    *destination = (int*) malloc(sizeof(int) * (*length));
    H5LTget_attribute_int(file_id, optionsObject, attributeName, *destination);

    pthread_mutex_unlock(&mutexHDF);
}

void dataproviderPrintInfo() {
    int maxwidth = 25;

    int numMultiStartRuns = getNumMultiStartRuns();
    logmessage(LOGLVL_INFO, "%*s: %d", maxwidth, "Num multistart optims", numMultiStartRuns);

    datapath path;

    for(int i = 0; i < numMultiStartRuns; ++i) {
        path.idxMultiStart = i;
        path.idxLocalOptimization = 0;
        int numStarts = getNumLocalOptimizationsForMultiStartRun(i);
        logmessage(LOGLVL_INFO, "%*s: %d", maxwidth, "Num starts", numStarts);
        logmessage(LOGLVL_INFO, "%*s: %d", maxwidth, "Genotypes", getNumGenotypes(path));
    }

    logmessage(LOGLVL_INFO, "%*s: %d", maxwidth, "Max iterations", getMaxIter());

}
