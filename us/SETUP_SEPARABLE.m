function Params = SETUP
% Taxing Human Capital
% @Author: Wenlan Luo
% Setup parameters

% Configure path
run '../SET_PATH.m';
set_dpopt;

% CMEX Compilation
% run CMEX file to preload memory
Pp = struct('form','MKLpp','breaks',{{[1 2 3 4]}},...
    'Values',[1 2 3 4],'coefs',[],'order',[4],...
    'Method',[],'ExtrapolationOrder',[],'thread',1,...
    'orient','curvefit');
WarmUpPp = myppual(Pp);

% Economic parameters
% These parameters can be overwriteen in subsequent routines, 
% or serve as initial guess of calibrating / optimization process

% earning shock
LogEpsilonMu = -0.029;
LogEpsilonSigma = 0.111;
EpsilonGridDispersion = 3;
EpsilonPts = 7;

% production function
Gamma = 0.36;
Delta = 0.0833;

% utility function
Beta = 0.9973;
Chi = 1.0217;
Sigma1 = 2;
Sigma2 = 3; % used only in separable utility function

% human capital elasticity and depreciation
Alpha = 0.5202;
Rho = 0;

% initial distribution of h and a
LogH1Mu = 0;
LogH1Sigma = 0.4019;
LogAMu = -1.6324;
LogASigma = 0.2713;
HARho = 0.7924;
APts = 7;
AGridDispersion = 3;  % max sigma dispersion in tauchen

% Demographic
Jw = 45;
Jr = 35;
J = Jw + Jr;
PopGrowthRate = 0.011;
% popg = 0;

% Death rate
SurvivalRate = csvread('alive.csv');
SurvivalRate = SurvivalRate(1:J);
% alive(:) = 1;

% Tax code
Tau0Gs = 0.258;
Tau1Gs = 0.768;
Tau2Gs = 2.6076;

% initial guess
Tau0K = 0.258;
Tau1K = 0.268;
Tau2K = 2.2543;
% Note(wenlan): Since I also found linear capital tax to be optimal,
% uncomment the following to shutdown capital income tax progressivity.
%{
Tau1K = 0;
Tau2K = 1;
%}

Tau0N = 0.258;
Tau1N = 0.768;
Tau2N = 2.2543;

% Note(wenlan): Result of total tax and factor tax should agree with each other, when
% tau2 = 0. Uncomment the following to test consistency
%{
Tau2Gs = 0;
Tau2K = 0;
Tau2N = 0;
%}

% social security and consumption tax
TauSs = 0.124;
TauC = 0.05;

% price
KLRatio = 4.7207;
Tss = 0.3976;
Tr = 0.0355;
r = Gamma * KLRatio^(Gamma-1) - Delta;
w = (1-Gamma) * KLRatio^Gamma;

EqX0 = [Beta Chi Tau2Gs Tss Tr];
BEqX0 = [Tau0N KLRatio Tss Tr];
BEqX0 = [0.3417    4.6599    0.3849    0.0323];
BEqLb = [0 1 0 0];
BEqUb = [1 9 Inf Inf];

GBench = 0.163174;

% Computation parameters
% Grids
% physical capital grid
KGrid = csvread('k_grid.csv');
KGrid = KGrid(1:6:end);
KMin = min(KGrid);
KPts = length(KGrid);
KMax= max(KGrid);

% human capital grid
HGrid = csvread('h_grid.csv') + 1e-1;
HGrid = HGrid(1:6:end);
HMin = min(HGrid);
HPts = length(HGrid);
HMax = max(HGrid);

% equiblrium parameters
TolEq = 1e-2;
TolOpt = 1e-6;
TolCon = 1e-6;
SpeedEqSmall = 0.1;
SpeedEqBig = 0.5;
IterMaxEq = 100;

% computation parameter
NumOfThreads = 8;
NSlots = getenv('NSLOTS');
if (strcmp(NSlots, '') == 0)
    NumOfThreads = str2double(NSlots)
end

% common used optimization boundary
KPrimeMin = KMin;
SMin = 1e-6;
NMin = 1e-6;
LMin = 1e-6;

CMin = 1e-6;

KPrimeMax = KMax;
SMax = 1 - NMin - LMin;
NMax = 1 - SMin - LMin;

% random number
NumOfAgents = 1e4;
Seed = 0823;

% calibrate bound
CaliX0 = [Alpha LogH1Sigma LogAMu LogASigma HARho];
CaliLb = [0.1 0 -2 0 -1];
CaliUb = [0.9 1 -1 0.3 1];
CaliScale = [0.1 0.1 0.1 0.01 0.1];

% tax search bound
TauX0 = [Tau0K Tau1K Tau2K Tau1N Tau2N];
TauX0 = [0.2467    0.0010    1.4014    3.0494   39.6992];
Tau0K = 0.2467;
Tau1K = 0.0010;
Tau2K = 1.4014;
Tau1N = 3.0494;
Tau2N = 39.6992;

TauLb = [0 1e-3 1 1e-3 1];
TauUb = [0.6 inf inf inf inf];
TauScale = [0.1 0.5 1 0.5 1];

% Assume linear tax on capital income
TauReducedX0 = [Tau0K Tau1N Tau2N];
TauReducedLb = [0 1e-3 0];
TauReducedUb = [0.7 inf inf];
TauScale = [0.1 0.5 1];

% pack local variables into structure
Params = v2struct();
end
