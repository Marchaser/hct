function tp = GS_PRIME(inc, tau0, tau1, tau2)
% compute marginal tax rate
bracket1 = 1 + tau2 * inc.^tau1;
bracket2 = 1 - bracket1.^(-1/tau1 - 1);
tp = tau0 * bracket2;
end