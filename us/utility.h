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
    PUSHN(a);
    PUSHN(k);
    PUSHN(h);
    PUSHN(kinc_before_tax);

#undef PUSH
#undef PUSHN
#undef PUSHARRAY
#undef PUSHNARRAY

#define MAX(a,b) (a>b) ? a : b
#define MIN(a,b) (a<b) ? a : b

    // controls
    double s = x[0];
    double n = x[1];
    double kp = x[2];

    // hp
    double hp_incr;
    double dhp_ds;
    if (s >= 1e-6) {
        hp_incr = a*pow(h*s, alpha);
        dhp_ds = alpha*hp_incr / s;
    }
    else {
        // fit a straight line connecting 0 and 1e-3
        double s_low = 1e-6;
        double hp_incr_at_low = a*pow(h*s_low, alpha);
        double slope = hp_incr_at_low / s_low;
        hp_incr = s*slope;
        dhp_ds = slope;
    }
    double hp = h*(1 - rho) + hp_incr;

    for (int i = 0; i < epsilon_pts; ++i) {
        stateFuture[i][0] = kp;
        stateFuture[i][1] = hp * epsilon_grid[i];
    }

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
    double ninc = (1 - tau_ss) * w*n*h;
    double y = ninc + kinc_before_tax;
    double tax;
    double dtax_dy;
    if (y >= 1e-6) {
        tax = tau0*(y - pow(pow(y, -tau1) + tau2, -1 / tau1));
        dtax_dy = tau0*(1 - pow(1 + tau2*pow(y, tau1), -1 / tau1 - 1));
    }
    else {
        // fit a straight line connecting 0 and 1e-3
        double y_low = 1e-6;
        double tax_low = tau0*(y_low - pow(pow(y_low, -tau1) + tau2, -1 / tau1));
        double slope = tax_low / y_low;
        tax = slope*y;
        dtax_dy = slope;
    }

    double wealth = k + Tr + y - tax;

    double c = (wealth - kp) / (1 + tau_c);
    double l = 1 - n - s;

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

    double dninc_dn = (1 - tau_ss) * (1 - dtax_dy) * w * h;

    if (f) {
        v = u;
        for (int i = 0; i < epsilon_pts; ++i) {
            v += vFuture[i] * epsilon_p[i];
        }
    }

    if (grad) {
        grad[0] = -dudl;
        grad[1] = dudc * dninc_dn / (1 + tau_c) - dudl;
        grad[2] = -dudc / (1 + tau_c);
        for (int i = 0; i < epsilon_pts; ++i) {
            grad[0] += gFuture[i][1] * dhp_ds*epsilon_grid[i] * epsilon_p[i];
            grad[2] += gFuture[i][0] * epsilon_p[i];
        }
    }

    if (cons) {
        cons[0] = c;
    }

    if (consgradRaw) {
        consgrad[0][0] = 0;
        consgrad[0][1] = dninc_dn / (1 + tau_c);
        consgrad[0][2] = -1 / (1 + tau_c);
    }

    USER_FUNC_RETURN;
}

