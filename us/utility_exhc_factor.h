/** define user function here:
 *
 * Please following the two steps to write user function:
 * 1. In first step, take x[] and data[] input, output futureState[][] which will be used in interpolation
 * 2. In second step, take x[], data[], vFuture[] and gFuture[] as input, output v, grad[]
 *  */
USER_FUNC_HEAD
{
    /**
     * preprocessor, do not edit
     */
    USER_FUNC_PRE;

    /**
     * step 1:
     * input: x[], data[]
     * output: stateFuture[][] of length [nvec, stateDim]
     */
    int i = 0;
#define PUSH(var) var = data[i++]
#define PUSHN(var) double var = data[i++]
#define PUSHARRAY(var, length) var = &data[i]; i+=length
#define PUSHNARRAY(var, length) double* var = &data[i]; i+=length
    int ii = 0;
#define GPUSH(var) var = g_shared_data[ii++]
#define GPUSHN(var) double var = g_shared_data[ii++]
#define GPUSHARRAY(var, length) var = &g_shared_data[ii]; ii+=length
#define GPUSHNARRAY(var, length) double* var = &g_shared_data[ii]; ii+=length

    GPUSHN(r);
    GPUSHN(w);
    GPUSHN(Tr);
    GPUSHN(tau0);
    GPUSHN(tau1);
    GPUSHN(tau2);
    GPUSHN(tau_ss);
    GPUSHN(tau_c);
    GPUSHN(alpha);
    GPUSHN(rho);
    GPUSHN(chi);
    GPUSHN(sigma1);
    GPUSHN(sigma2);
    GPUSHN(d_epsilon_pts);
    int epsilon_pts = (int)d_epsilon_pts;
    GPUSHNARRAY(epsilon_grid, epsilon_pts);
    GPUSHNARRAY(epsilon_p, epsilon_pts);


    // states
    PUSHN(k);
    PUSHN(sbar);
    PUSHN(hbar);
    PUSHN(kinc_after_tax);

#undef PUSH
#undef PUSHN
#undef PUSHARRAY
#undef PUSHNARRAY

#define MAX(a,b) (a>b) ? a : b
#define MIN(a,b) (a<b) ? a : b

    // controls
    double n = x[0];
    double kp = x[1];

    stateFuture[0][0] = kp;

    /**
     * interpolate future state, do not edit
     */
    USER_FUNC_INTERP;

    /**
     * step 2:
     * input: x[], data[], vFuture[], gFuture[][]
     * output:
     *    v: scalar value eavaluated at x
     *    grad[] of size controlDim: gradient of v evaluated at x
     *    cons[] of size nonlin: evaluation of constraint
     *    consgrad[][] of size [nonlin, controlDim]: jacobian of cons
     */

    // tax
    double ninc = (1 - tau_ss) * w*n*hbar;
    double ntax;
    double dntax_dninc;
    if (ninc >= 1e-6) {
        ntax = tau0*(ninc - pow(pow(ninc, -tau1) + tau2, -1 / tau1));
        dntax_dninc = tau0*(1 - pow(1 + tau2*pow(ninc, tau1), -1 / tau1 - 1));
    }
    else {
        // fit a straight line connecting 0 and 1e-3
        double ninc_low = 1e-6;
        double tax_low = tau0*(ninc - pow(pow(ninc_low, -tau1) + tau2, -1 / tau1));
        double slope = tax_low / ninc_low;
        ntax = slope*ninc;
        dntax_dninc = slope;
    }

    double wealth = k + Tr + kinc_after_tax + ninc - ntax;

    double c = (wealth - kp) / (1 + tau_c);
    double l = 1 - n - sbar;

    double inf = 1e20;

    // u
#ifdef NON_SEPARABLE
    double cl_comp = (c > 0 && l > 0) ? pow(c, chi) * pow(l, 1 - chi) : -inf;
    double pow_part = pow(cl_comp, 1 - sigma1);
    double u = (c > 0 && l > 0) ? pow_part / (1 - sigma1) : -inf;
    double du_dcl_comp = pow_part / cl_comp;
    double dudc = (c > 0) ? du_dcl_comp * chi * cl_comp / c : inf;
    double dudl = (l > 0) ? du_dcl_comp * (1 - chi) * cl_comp / l : inf;
#else
    double u = (c > 0 && l > 0) ? pow(c, 1 - sigma1) / (1 - sigma1) + chi*pow(l, 1 - sigma2) / (1 - sigma2) : -inf;
    double dudc = (c > 0) ? pow(c, -sigma1) : inf;
    double dudl = (l > 0) ? chi*pow(l, -sigma2) : inf;
#endif

    double dninc_dn = (1 - tau_ss) * (1 - dntax_dninc) * w * hbar;

    if (f) {
        v = u + vFuture[0];
    }

    if (grad) {
        grad[0] = dudc * dninc_dn / (1 + tau_c) - dudl;
        grad[1] = -dudc / (1 + tau_c) + gFuture[0][0];
    }

    if (cons) {
        cons[0] = c;
    }

    if (consgradRaw) {
        consgrad[0][0] = dninc_dn / (1 + tau_c);
        consgrad[0][1] = -1 / (1 + tau_c);
    }

    USER_FUNC_RETURN;
}

