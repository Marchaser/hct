function NewParams = OVERWRITE(varargin)
Params = varargin{1};
v2struct(Params);

if nargin>1
    Profile = varargin{2};
    switch(Profile)
        case 'Benchmark'
            BenchResultFileName = 'BenchResult';
            GBench = 0.1440;
            CaliResult = [0.3747    0.4342   -1.6462    0.2682    0.9189];
            CaliResultEqX = [1.0128    0.3209    2.7342    0.3506    0.0296];
        case 'NoBc'
            BenchResultFileName = 'NoBc';
            KGrid = KGrid-1; % extend to the left
            KPts = length(KGrid);
            KMax = max(KGrid);
            KMin = min(KGrid);
            KPrimeMin = KMin; % can borrow.
            GBench = 0.144413;
            CaliResult = [0.37051      0.43578      -1.6604      0.26538      0.94194];
            CaliResultEqX = [1.0158      0.32039       2.6525      0.35173     0.027947];
        case 'NoIdio'
%             LogEpsilonMu = LogEpsilonMu + 0.5*LogEpsilonSigma^2;
%             LogEpsilonSigma = 1e-12;
%             CaliResult = [0.3747    0.4342   -1.6462    0.2682    0.9189];
%             CaliResultEqX = [1.0128    0.3209    2.7342    0.3506    0.0296];
%             BenchResultFileName = 'NoIdio';
            LogEpsilonSigma = 1e-12;
            LogEpsilonMu = -0.025713;
            BenchResultFileName = 'NoIdio';
            GBench = 0.134697;
            CaliResult = [0.3747    0.4342   -1.6462    0.2682    0.9189];
            CaliResultEqX = [1.0214      0.30506       3.0404      0.32802     0.029508];
        case 'NoType'
            LogH1Mu = -0.040224;
            BenchResultFileName = 'NoType';
            GBench = 0.122556;
            CaliResult = [0.3747    1e-12   -1.6723    1e-12    0.9189];
            CaliResultEqX = [1.0136      0.31359       3.2928      0.29839     0.024524];
    end
end

if nargin>2
    OptimTau = varargin{3};
    switch(OptimTau)
        case 'Base'
            TauX = [0.45331       0 0       11.476    4.674e+08];
            BEqX = [0.17825        4.156      0.34902     0.027595];
            TauX0 = [Tau0K Tau1K Tau2K Tau1N Tau2N];
            BEqX0 = [Tau0N KLRatio Tss Tr];
        case 'Exhc'
            TauX = [0.36334     0   0       10.832   4.4328e+08];
            BEqX = [0.22246 4.4245      0.34351     0.02754];
            TauX0 = [Tau0K Tau1K Tau2K Tau1N Tau2N];
            BEqX0 = [Tau0N KLRatio Tss Tr];
        case 'Exn'
%             TauX = [0.3317     0   0    11.35   5.0589e+08];
%             BEqX = [0.34383       4.8756      0.30559     0.02366];
            GBench = 0.130853;
            CaliResult = [0.484      0.43836      -1.4546       0.1938      0.89502];
            CaliResultEqX = [1.0158      0.25078       3.0641      0.31851     0.026823];
            NMin = 0.2641;
            NMax = 0.2641;
            TauX0 = [Tau0K Tau1K Tau2K Tau1N Tau2N];
            BEqX0 = [Tau0N KLRatio Tss Tr];
%             KShrinkFactor = 5;
%             HShrinkFactor = 5;
%             TauX0 = [0.3317     0   0    11.35   5.0589e+08];
%             BEqX0 = [0.34383       4.8756      0.30559     0.02366];
        case 'ExAll'
            TauX = [-0.84803     0.057516       164.16       11.625       9610.4];
            BEqX = [1 6.345      0.35433     0.031596];
            
            BenchResultFileName = 'Exn';
            GBench = 0.130853;
            CaliResult = [0.484      0.43836      -1.4546       0.1938      0.89502];
            CaliResultEqX = [1.0158      0.25078       3.0641      0.31851     0.026823];
            NMin = 0.2641;
            NMax = 0.2641;
            TauX0 = [Tau0K Tau1K Tau2K Tau1N Tau2N];
            BEqX0 = [Tau0N KLRatio Tss Tr];
    end
end

if isequal(exist('CaliResult','var'),1)
    Alpha = CaliResult(1);
    LogH1Sigma = CaliResult(2);
    LogAMu = CaliResult(3);
    LogASigma = CaliResult(4);
    HARho = CaliResult(5);
    CaliX0 = CaliResult;
end

if isequal(exist('CaliResultEqX','var'),1)
    Beta = CaliResultEqX(1);
    Chi = CaliResultEqX(2);
    Tau2Gs = CaliResultEqX(3);
    Tss = CaliResultEqX(4);
    Tr = CaliResultEqX(5);
    EqX0 = CaliResultEqX;
end

if isequal(exist('TauX','var'),1)
    Tau0K = TauX(1);
    Tau1K = TauX(2);
    Tau2K = TauX(3);
    Tau1N = TauX(4);
    Tau2N = TauX(5);
end

if isequal(exist('BEqX','var'),1)
    Tau0N = BEqX(1);
    KLRatio = BEqX(2);
    Tss = BEqX(3);
    Tr = BEqX(4);
    r = Gamma * KLRatio^(Gamma-1) - Delta;
    w = (1-Gamma) * KLRatio^Gamma;
end

clear Params;
NewParams = v2struct;
end