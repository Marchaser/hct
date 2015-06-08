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

    GPUSHN(tau_c);
    GPUSHN(sigma1);
    GPUSHN(chi);

    PUSHN(budget_after_tax);

#undef PUSH
#undef PUSHN
#undef PUSHARRAY
#undef PUSHNARRAY

#define MAX(a,b) (a>b) ? a : b
#define MIN(a,b) (a<b) ? a : b

    // controls
    double kp = x[0];
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

    double y = budget_after_tax;

    double c = (y - kp) / (1 + tau_c);

    if (f) {
#ifdef NON_SEPARABLE
        v = (c > 0) ? pow(c, chi*(1 - sigma1)) / (1 - sigma1) + vFuture[0] : -1e8;
#else
        v = (c > 0) ? pow(c, 1 - sigma1) / (1 - sigma1) + vFuture[0] : -1e8;
#endif
    }

    if (grad) {
    }

    if (cons) {
    }

    if (consgradRaw) {
    }

    USER_FUNC_RETURN;
}

