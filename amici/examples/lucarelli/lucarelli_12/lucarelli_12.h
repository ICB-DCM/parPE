#ifndef _amici_TPL_MODELNAME_h
#define _amici_TPL_MODELNAME_h
#include <cmath>
#include <memory>
#include "amici/defines.h"
#include <sundials/sundials_sparse.h> //SlsMat definition
#include "amici/solver_cvodes.h"
#include "amici/model_ode.h"

namespace amici {
class Solver;
}

/**
 * @brief Wrapper function to instantiate the linked Amici model without knowing the name at compile time.
 * @return
 */
extern void J_lucarelli_12(realtype *J, const realtype t, const realtype *x, const double *p, const double *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void JB_lucarelli_12(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
extern void JDiag_lucarelli_12(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void JSparse_lucarelli_12(SlsMat JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx);
extern void JSparseB_lucarelli_12(SlsMat JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
extern void Jv_lucarelli_12(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *v, const realtype *w, const realtype *dwdx);
extern void JvB_lucarelli_12(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *vB, const realtype *w, const realtype *dwdx);
extern void Jy_lucarelli_12(double *nllh, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my);
extern void dJydsigma_lucarelli_12(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my);
extern void dJydy_lucarelli_12(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my);
extern void dwdp_lucarelli_12(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dwdx_lucarelli_12(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void dxdotdp_lucarelli_12(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp);
extern void dydx_lucarelli_12(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void dydp_lucarelli_12(double *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip);
extern void dsigmaydp_lucarelli_12(double *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip);
extern void qBdot_lucarelli_12(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp);
extern void sigmay_lucarelli_12(double *sigmay, const realtype t, const realtype *p, const realtype *k);
extern void sxdot_lucarelli_12(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *sx, const realtype *w, const realtype *dwdx, const realtype *J, const realtype *dxdotdp);
extern void w_lucarelli_12(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
extern void x0_lucarelli_12(realtype *x0, const realtype t, const realtype *p, const realtype *k);
extern void xBdot_lucarelli_12(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx);
extern void xdot_lucarelli_12(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w);
extern void y_lucarelli_12(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);

/**
 * @brief AMICI-generated model subclass.
 */
class Model_lucarelli_12 : public amici::Model_ODE {
public:
    /**
     * @brief Default constructor.
     */
    Model_lucarelli_12()
        : amici::Model_ODE(
              33, // nx
              33, // nxtrue
              12, // ny
              12, // nytrue
              0, // nz
              0, // nztrue
              0, // nevent
              1, // nobjective
              57, // nw
              93, // ndwddx
              129, // ndwdp
              141, // nnz
              33, // ubw
              33, // lbw
              amici::AMICI_O2MODE_NONE, // o2mode
              std::vector<realtype>{0.00598343088790046, 0.286397320264704, 0.0489949146416608, 0.379720142751521, 0.0141389661164018, 0.000797280887701997, 0.0, 0.0, 0.0462523981984203, 0.026381058989095, 0.00455788978022206, 0.461699170803527, 0.073852775177442, 0.0304255473876269, 0.0, 0.747200749545747, 0.407329851468312, 0.0112459404708192, 0.00714636599139809, 0.0, 0.0, 0.0, 0.018441582221635, 0.035379693264711, 0.00378630404521753, 0.0111673720417551, 0.000268942132758932, 0.0, 0.0, 0.0868803385899104, 0.101419105775767, 0.000802160834341515, 0.0, 1.03661055568903, 6.10224669858575, 8.20035989079744, 1.4392548717489, 9.67035015330535, 0.124380538819967, 998.299919973209, 0.135655043663672, 20.6442242166861, 3.63157061548523, 0.0, 0.0, 36.6724104315148, 626.768522107231, 55.9534537604082, 0.0, 0.0, 0.0175854489108708, 0.000815621975867545, 295.971427530835, 0.0, 89.9866382703511, 999.921986303898, 999.998960700206, 19.0317907364409, 218.356690800299, 0.992827361229694, 9.30152789301652, 0.0, 0.0, 0.0, 320.763519415643, 999.981624754706, 0.00328197841814657, 0.0, 0.0, 0.0, 0.47739150796461, 0.00985305306597475, 0.00851371884219428, 0.0140190729281921, 0.0125133830064198, 0.000284224503418874, 0.0, 0.0, 1.36603583936842, 2.43594567888689, 0.00359915570352007, 0.0, 0.000474849473374749, 0.0, 1.29976503414042, 0.062323573525662, 0.13414856027346, 0.0149504509077553, 1.83614829885072, 0.0, 0.0, 0.148655616115831, 0.000403100938582095, 8.56307475687022e-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0404498417896043, 0.140066630704836, 0.133558955086407, 0.0961378774278948, 0.104874995503069, 0.163608622879931, 0.116787493223069, 0.0733384767144999, 0.101509175364595, 0.173716659133795, 0.0699132545087422, 0.125435007986767, 0.0766684551714607}, // dynamic parameters
              std::vector<realtype>{1.0, 142.777172927249, 16.258584970194, 67.0508150691882}, // fixedParameters
              std::vector<int>{}, // plist
              std::vector<realtype>(33,0.0), // idlist
              std::vector<int>{} // z2event
              )
    {}

    /**
     * @brief Clone this model instance.
     * @return A deep copy of this instance.
     */
    virtual amici::Model* clone() const override { return new Model_lucarelli_12(*this); }

    /** model specific implementation for fJ
     * @param J Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJ(realtype *J, const realtype t, const realtype *x, const double *p, const double *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        J_lucarelli_12(J, t, x, p, k, h, w, dwdx);
    }

    /** model specific implementation for fJB
     * @param JB Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param xB Vector with the adjoint states
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJB(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) override {
        JB_lucarelli_12(JB, t, x, p, k, h, xB, w, dwdx);
    }

    /** model specific implementation for fJDiag
     * @param JDiag Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJDiag(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        JDiag_lucarelli_12(JDiag, t, x, p, k, h, w, dwdx);
    }

    /** model specific implementation for fJSparse
     * @param JSparse Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJSparse(SlsMat JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) override {
        JSparse_lucarelli_12(JSparse, t, x, p, k, h, w, dwdx);
    }

    /** model specific implementation for fJSparseB
     * @param JSparseB Matrix to which the Jacobian will be written
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param xB Vector with the adjoint states
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJSparseB(SlsMat JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) override {
        JSparseB_lucarelli_12(JSparseB, t, x, p, k, h, xB, w, dwdx);
    }

    /** model specific implementation of fJrz
     * @param nllh regularization for event measurements z
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     **/
    virtual void fJrz(double *nllh, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) override {
    }

    /** model specific implementation for fJv
     * @param Jv Matrix vector product of J with a vector v
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param v Vector with which the Jacobian is multiplied
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJv(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *v, const realtype *w, const realtype *dwdx) override {
        Jv_lucarelli_12(Jv, t, x, p, k, h, v, w, dwdx);
    }

    /** model specific implementation for fJvB
     * @param JvB Matrix vector product of JB with a vector v
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param xB Vector with the adjoint states
     * @param vB Vector with which the Jacobian is multiplied
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fJvB(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *vB, const realtype *w, const realtype *dwdx) override {
        JvB_lucarelli_12(JvB, t, x, p, k, h, xB, vB, w, dwdx);
    }

    /** model specific implementation of fJy
     * @param nllh negative log-likelihood for measurements y
     * @param iy output index
     * @param p parameter vector
     * @param k constant vector
     * @param y model output at timepoint
     * @param sigmay measurement standard deviation at timepoint
     * @param my measurements at timepoint
     **/
    virtual void fJy(double *nllh, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) override {
        Jy_lucarelli_12(nllh, iy, p, k, y, sigmay, my);
    }

    /** model specific implementation of fJz
     * @param nllh negative log-likelihood for event measurements z
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurements at timepoint
     **/
    virtual void fJz(double *nllh, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) override {
    }

    /** model specific implementation of fdJrzdsigma
     * @param dJrzdsigma Sensitivity of event penalization Jrz w.r.t.
     * standard deviation sigmaz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param rz model root output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     **/
    virtual void fdJrzdsigma(double *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) override {
    }

    /** model specific implementation of fdJrzdz
     * @param dJrzdz partial derivative of event penalization Jrz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param rz model root output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     **/
    virtual void fdJrzdz(double *dJrzdz, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) override {
    }

    /** model specific implementation of fdJydsigma
     * @param dJydsigma Sensitivity of time-resolved measurement
     * negative log-likelihood Jy w.r.t. standard deviation sigmay
     * @param iy output index
     * @param p parameter vector
     * @param k constant vector
     * @param y model output at timepoint
     * @param sigmay measurement standard deviation at timepoint
     * @param my measurement at timepoint
     **/
    virtual void fdJydsigma(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) override {
        dJydsigma_lucarelli_12(dJydsigma, iy, p, k, y, sigmay, my);
    }

    /** model specific implementation of fdJydy
     * @param dJydy partial derivative of time-resolved measurement negative log-likelihood Jy
     * @param iy output index
     * @param p parameter vector
     * @param k constant vector
     * @param y model output at timepoint
     * @param sigmay measurement standard deviation at timepoint
     * @param my measurement at timepoint
     **/
    virtual void fdJydy(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) override {
        dJydy_lucarelli_12(dJydy, iy, p, k, y, sigmay, my);
    }

    /** model specific implementation of fdJzdsigma
     * @param dJzdsigma Sensitivity of event measurement
     * negative log-likelihood Jz w.r.t. standard deviation sigmaz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurement at timepoint
     **/
    virtual void fdJzdsigma(double *dJzdsigma, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) override {
    }

    /** model specific implementation of fdJzdz
     * @param dJzdz partial derivative of event measurement negative log-likelihood Jz
     * @param iz event output index
     * @param p parameter vector
     * @param k constant vector
     * @param z model event output at timepoint
     * @param sigmaz event measurement standard deviation at timepoint
     * @param mz event measurement at timepoint
     **/
    virtual void fdJzdz(double *dJzdz, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) override {
    }

    /** model specific implementation of fdeltasx
     * @param deltaqB sensitivity update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ip sensitivity index
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param xB adjoint state
     **/
    virtual void fdeltaqB(double *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB) override {
    }

    /** model specific implementation of fdeltasx
     * @param deltasx sensitivity update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param w repeating elements vector
     * @param ip sensitivity index
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param sx state sensitivity
     * @param stau event-time sensitivity
     **/
    virtual void fdeltasx(double *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau) override {
    }

    /** model specific implementation of fdeltax
     * @param deltax state update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     **/
    virtual void fdeltax(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old) override {
    }

    /** model specific implementation of fdeltaxB
     * @param deltaxB adjoint state update
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ie event index
     * @param xdot new model right hand side
     * @param xdot_old previous model right hand side
     * @param xB current adjoint state
     **/
    virtual void fdeltaxB(double *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB) override {
    }

    /** model specific implementation of fdrzdp
     * @param drzdp partial derivative of root output rz w.r.t. model parameters p
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     **/
    virtual void fdrzdp(double *drzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {
    }

    /** model specific implementation of fdrzdx
     * @param drzdx partial derivative of root output rz w.r.t. model states x
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void fdrzdx(double *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
    }

    /** model specific implementation of fsigmay
     * @param dsigmaydp partial derivative of standard deviation of measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fdsigmaydp(double *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip) override {
        dsigmaydp_lucarelli_12(dsigmaydp, t, p, k, ip);

    }

    /** model specific implementation of fsigmaz
     * @param dsigmazdp partial derivative of standard deviation of event measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fdsigmazdp(double *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip) override {
    }

    /** model specific implementation of dwdp
     * @param dwdp Recurring terms in xdot, parameter derivative
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     */
    virtual void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dwdp_lucarelli_12(dwdp, t, x, p, k, h, w);
    }

    /** model specific implementation of dwdx
     * @param dwdx Recurring terms in xdot, state derivative
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     */
    virtual void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        dwdx_lucarelli_12(dwdx, t, x, p, k, h, w);
    }

    /** model specific implementation of fdxdotdp
     * @param dxdotdp partial derivative xdot wrt p
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param ip parameter index
     * @param w vector with helper variables
     * @param dwdp derivative of w wrt p
     */
    virtual void fdxdotdp(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp) override {
        dxdotdp_lucarelli_12(dxdotdp, t, x, p, k, h, ip, w, dwdp);
    }

    /** model specific implementation of fdydx
     * @param dydx partial derivative of observables y w.r.t. model states x
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void fdydx(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        dydx_lucarelli_12(dydx, t, x, p, k, h);
    }

    /** model specific implementation of fdydp
     * @param dydp partial derivative of observables y w.r.t. model parameters p
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     **/
    virtual void fdydp(double *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {
        dydp_lucarelli_12(dydp, t, x, p, k, h, ip);
    }

    /** model specific implementation of fdzdp
     * @param dzdp partial derivative of event-resolved output z w.r.t. model parameters p
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param ip parameter index w.r.t. which the derivative is requested
     **/
    virtual void fdzdp(double *dzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) override {
    }

    /** model specific implementation of fdzdx
     * @param dzdx partial derivative of event-resolved output z w.r.t. model states x
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void fdzdx(double *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
    }

    /** model specific implementation for fqBdot
     * @param qBdot adjoint quadrature equation
     * @param ip sensitivity index
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param xB Vector with the adjoint states
     * @param w vector with helper variables
     * @param dwdp derivative of w wrt p
     **/
    virtual void fqBdot(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp) override {
        qBdot_lucarelli_12(qBdot, ip, t, x, p, k, h, xB, w, dwdp);
    }

    /** model specific implementation for froot
     * @param root values of the trigger function
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     **/
    virtual void froot(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
    }

    /** model specific implementation of frz
     * @param rz value of root function at current timepoint (non-output events not included)
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void frz(double *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
    }

    /** model specific implementation of fsigmay
     * @param sigmay standard deviation of measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fsigmay(double *sigmay, const realtype t, const realtype *p, const realtype *k) override {
        sigmay_lucarelli_12(sigmay, t, p, k);
    }

    /** model specific implementation of fsigmaz
     * @param sigmaz standard deviation of event measurements
     * @param t current time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fsigmaz(double *sigmaz, const realtype t, const realtype *p, const realtype *k) override {
    }

    /** model specific implementation of fsrz
     * @param srz Sensitivity of rz, total derivative
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param sx current state sensitivity
     * @param h heavyside vector
     * @param ip sensitivity index
     **/
    virtual void fsrz(double *srz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip) override {
    }

    /** model specific implementation of fstau
     * @param stau total derivative of event timepoint
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param sx current state sensitivity
     * @param ip sensitivity index
     * @param ie event index
     **/
    virtual void fstau(double *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip, const int ie) override {
    }

    /** model specific implementation of fsx0
     * @param sx0 initial state sensitivities
     * @param t initial time
     * @param x0 initial state
     * @param p parameter vector
     * @param k constant vector
     * @param ip sensitivity index
     **/
    virtual void fsx0(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip) override {
    }

    /** model specific implementation of fsxdot
     * @param sxdot sensitivity rhs
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param ip parameter index
     * @param sx Vector with the state sensitivities
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     * @param J jacobian
     * @param dxdotdp parameter derivative of residual function
     */
    virtual void fsxdot(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *sx, const realtype *w, const realtype *dwdx, const realtype *J, const realtype *dxdotdp) override {
        sxdot_lucarelli_12(sxdot, t, x, p, k, h, ip, sx, w, dwdx, J, dxdotdp);
    }

    /** model specific implementation of fsz
     * @param sz Sensitivity of rz, total derivative
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     * @param sx current state sensitivity
     * @param ip sensitivity index
     **/
    virtual void fsz(double *sz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip) override {
    }

    /** model specific implementation of fw
     * @param w Recurring terms in xdot
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     */
    virtual void fw(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        w_lucarelli_12(w, t, x, p, k, h);
    }

    /** model specific implementation of fx0
     * @param x0 initial state
     * @param t initial time
     * @param p parameter vector
     * @param k constant vector
     **/
    virtual void fx0(realtype *x0, const realtype t, const realtype *p, const realtype *k) override {
        x0_lucarelli_12(x0, t, p, k);
    }

    /** model specific implementation for fxBdot
     * @param xBdot adjoint residual function
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param xB Vector with the adjoint states
     * @param w vector with helper variables
     * @param dwdx derivative of w wrt x
     **/
    virtual void fxBdot(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) override {
        xBdot_lucarelli_12(xBdot, t, x, p, k, h, xB, w, dwdx);
    }

    /** model specific implementation for fxdot
     * @param xdot residual function
     * @param t timepoint
     * @param x Vector with the states
     * @param p parameter vector
     * @param k constants vector
     * @param h heavyside vector
     * @param w vector with helper variables
     **/
    virtual void fxdot(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) override {
        xdot_lucarelli_12(xdot, t, x, p, k, h, w);
    }

    /** model specific implementation of fy
     * @param y model output at current timepoint
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void fy(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
        y_lucarelli_12(y, t, x, p, k, h);
    }

    /** model specific implementation of fz
     * @param z value of event output
     * @param ie event index
     * @param t current time
     * @param x current state
     * @param p parameter vector
     * @param k constant vector
     * @param h heavyside vector
     **/
    virtual void fz(double *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) override {
    }



    /**
     * @brief Get names of the model parameters
     * @return the names
     */
    virtual std::vector<std::string> getParameterNames() const { return std::vector<std::string> {"Rec_act",
"S_dephos",
"S_dephosphos",
"S_phos",
"geneA_act1",
"geneA_act2",
"geneA_act3",
"geneA_inh1",
"geneA_inh2",
"geneA_inh3",
"geneA_turn",
"geneB_act1",
"geneB_act2",
"geneB_act3",
"geneB_inh1",
"geneB_inh2",
"geneB_inh3",
"geneB_turn",
"geneC_act1",
"geneC_act2",
"geneC_act3",
"geneC_inh1",
"geneC_inh2",
"geneC_inh3",
"geneC_turn",
"geneD_act1",
"geneD_act2",
"geneD_act3",
"geneD_inh1",
"geneD_inh2",
"geneD_inh3",
"geneD_turn",
"geneE_act1",
"geneE_act2",
"geneE_act3",
"geneE_inh1",
"geneE_inh2",
"geneE_inh3",
"geneE_turn",
"geneF_act1",
"geneF_act2",
"geneF_act3",
"geneF_inh1",
"geneF_inh2",
"geneF_inh3",
"geneF_turn",
"geneG_act1",
"geneG_act2",
"geneG_act3",
"geneG_inh1",
"geneG_inh2",
"geneG_inh3",
"geneG_turn",
"geneH_act1",
"geneH_act2",
"geneH_act3",
"geneH_inh1",
"geneH_inh2",
"geneH_inh3",
"geneH_turn",
"geneI_act1",
"geneI_act2",
"geneI_act3",
"geneI_inh1",
"geneI_inh2",
"geneI_inh3",
"geneI_turn",
"geneJ_act1",
"geneJ_act2",
"geneJ_act3",
"geneJ_inh1",
"geneJ_inh2",
"geneJ_inh3",
"geneJ_turn",
"geneK_act1",
"geneK_act2",
"geneK_act3",
"geneK_inh1",
"geneK_inh2",
"geneK_inh3",
"geneK_turn",
"geneL_act1",
"geneL_act2",
"geneL_act3",
"geneL_inh1",
"geneL_inh2",
"geneL_inh3",
"geneL_turn",
"init_Rec",
"k_223",
"k_224",
"k_233",
"k_234",
"k_244",
"k_334",
"k_344",
"k_on_u",
"kdiss_SS",
"khomo2",
"khomo3",
"khomo4",
"pRec_degind",
"sd_Bmp4_nExpID100",
"sd_Cxcl15_nExpID100",
"sd_Dnmt3a_nExpID100",
"sd_Dusp5_nExpID100",
"sd_Jun_nExpID100",
"sd_Klf10_nExpID100",
"sd_Pdk4_nExpID100",
"sd_Ski_nExpID100",
"sd_Skil_nExpID100",
"sd_Smad7_nExpID100",
"sd_Sox4_nExpID100",
"sd_Tgfa_nExpID100",}; }

    /**
     * @brief Get names of the model states
     * @return the names
     */
    virtual std::vector<std::string> getStateNames() const { return std::vector<std::string> {"TGFb",
"Rec",
"TGFb_pRec",
"S2",
"S3",
"S4",
"S2_S4_S4",
"ppS2_ppS2_ppS2",
"ppS3_ppS3_ppS3",
"S4_S4_S4",
"pS2",
"pS3",
"ppS2",
"ppS3",
"ppS2_ppS2_S4",
"ppS2_ppS2_ppS3",
"ppS2_ppS3_ppS3",
"ppS3_ppS3_S4",
"ppS2_ppS3_S4",
"ppS3_S4_S4",
"ppS2_S4_S4",
"geneA",
"geneB",
"geneC",
"geneD",
"geneE",
"geneF",
"geneG",
"geneH",
"geneI",
"geneJ",
"geneK",
"geneL",}; }

    /**
     * @brief Get names of the fixed model parameters
     * @return the names
     */
    virtual std::vector<std::string> getFixedParameterNames() const { return std::vector<std::string> {"init_TGFb",
"init_S2",
"init_S3",
"init_S4",}; }

    /**
     * @brief Get names of the observables
     * @return the names
     */
    virtual std::vector<std::string> getObservableNames() const { return std::vector<std::string> {"observable_Ski",
"observable_Skil",
"observable_Dnmt3a",
"observable_Sox4",
"observable_Jun",
"observable_Smad7",
"observable_Klf10",
"observable_Bmp4",
"observable_Cxcl15",
"observable_Dusp5",
"observable_Tgfa",
"observable_Pdk4",}; }

};

#endif /* _amici_TPL_MODELNAME_h */
