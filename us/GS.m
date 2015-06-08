function tax = GS(y, tau0, tau1, tau2)
if tau1<=0
    tax = tau0*y;
    return;
end
bracket1 = y.^(-tau1) + tau2;
bracket2 = y - bracket1.^(-1/tau1);
tax = tau0 * bracket2;
% tax(y<=0) = tau0*y(y<=0);
tax(y<=0) = 0;
end