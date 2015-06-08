function [SimulateResult, ExitFlag] = SIMULATE(VfiResult, Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Simulate economy, return distribution

% unpack structure
v2struct(Params);
v2struct(VfiResult);

% initialize storage
% Note(wenlan): variable starting with lower case is individual variable,
% to be distinguished from the aggregate variable starting with upper case
SimulateResult.c = zeros(NumOfAgents, J);
SimulateResult.k = zeros(NumOfAgents, J+1);
SimulateResult.h = zeros(NumOfAgents, J);
SimulateResult.s = zeros(NumOfAgents, J);
SimulateResult.n = zeros(NumOfAgents, J);

SimulateResult.v = zeros(NumOfAgents, J);

SimulateResult.kinc = zeros(NumOfAgents, J);
SimulateResult.ninc = zeros(NumOfAgents, J);
SimulateResult.tax = zeros(NumOfAgents, J);
SimulateResult.tss = zeros(NumOfAgents, J);
SimulateResult.survival = ones(NumOfAgents, J);
SimulateResult.population = ones(NumOfAgents, J+1);

% initial distribution
SimulateResult.k(:,1) = K1;
SimulateResult.h(:,1) = H1;

% simulate young
% NOTE(wenlan): Construct value function spline. Keep in mind, vec has the memory order (age, epsilon, a)
%{
ValuePp = struct('form','pp','breaks',{{KGrid,HGrid}},...
    'Values',[],'coefs',[],'order',[4 4],...
    'Method',[],'ExtrapolationOrder',[],'thread',NumOfThreads,...
    'orient','curvefit');
ValuePp.Values = reshape(permute(ValueYoung, [4 1 2 3]), [], KPts, HPts);
ValuePp = myppual(ValuePp);
%}
ValuePp = tensor_pchip({KGrid,HGrid}, reshape(permute(ValueYoung, [4 1 2 3]), [], KPts, HPts));
ProblemYoung.pp = myppual(ValuePp);
ValueVecSize = [Jw+1 APts];

ProblemYoung.options.Algorithm = 'donlp2';
ProblemYoung.options.NumThreads = NumOfThreads;
ProblemYoung.options.TolX = TolOpt;
ProblemYoung.lb = [
    SMin * ones(1, NumOfAgents)
    NMin * ones(1, NumOfAgents)
    KPrimeMin * ones(1, NumOfAgents)
    ];

ProblemYoung.Aineq = repmat([1;1;0], 1, NumOfAgents);
ProblemYoung.blb = zeros(1, NumOfAgents);
ProblemYoung.bub = ones(1, NumOfAgents);
ProblemYoung.clb = CMin * ones(1, NumOfAgents);
ProblemYoung.cub = 1e20 * ones(1, NumOfAgents);

ProblemYoung.shared_data = [
    r; w; Tr;
    Tau0Gs; Tau1Gs; Tau2Gs;
    TauSs; TauC;
    Alpha; Rho;
    Chi; Sigma1; Sigma2;
    EpsilonPts;
    EpsilonGrid(:);
    EpsilonP(:);
    ];

X = [
    SMin * ones(1, NumOfAgents);
    NMin * ones(1, NumOfAgents);
    KPrimeMin * ones(1, NumOfAgents)
    ];

for j=1:1:Jw
  % compute future value interpolation vector function index
  VFutureIdx = repmat(sub2ind(ValueVecSize, (j+1)*ones(1, NumOfAgents), AIdx(:)'), EpsilonPts, 1);

  % capital income
  KInc = r * (SimulateResult.k(:,j) + Tr);
  KWealth = (1+r) * (SimulateResult.k(:,j) + Tr);
  
  % labor income when working full time
  MaxNInc = (1-TauSs) * NMax * w * SimulateResult.h(:,j);
  
  % maximum income after tax
  MaxBudget = KWealth + MaxNInc - GS(KInc + MaxNInc, Tau0Gs, Tau1Gs, Tau2Gs);

  % use this period distribution as data
  ProblemYoung.Data = [...
    AShock(:)';
    SimulateResult.k(:,j)';
    SimulateResult.h(:,j)';
    KInc(:)';
    ];

  ProblemYoung.ub = [
      SMax * ones(1, NumOfAgents)
	  NMax * ones(1, NumOfAgents)
	  min(MaxBudget(:)'-CMin, KPrimeMax)
      ];
  ProblemYoung.idx = VFutureIdx;
  
  ProblemYoung.x0 = X;
  [X, V, OptExitFlag] = dpopt_utility(ProblemYoung);

  % extract policy
  SimulateResult.s(:,j) = X(1,:);
  SimulateResult.n(:,j) = X(2,:);
  SimulateResult.k(:,j+1) = X(3,:);
  SimulateResult.v(:,j) = V;
  
  % next period h
  SimulateResult.h(:,j+1) = EpsilonShock(:,j) .* ...
      ((1-Rho)*SimulateResult.h(:,j) + AShock(:) .* (SimulateResult.s(:,j).*SimulateResult.h(:,j)).^Alpha);
  
  
  NInc = w * (1-TauSs) * SimulateResult.n(:,j) .* SimulateResult.h(:,j);
  SimulateResult.tss(:,j) = w*TauSs*SimulateResult.n(:,j) .* SimulateResult.h(:,j);
  SimulateResult.tax(:,j) = GS(KInc + NInc, Tau0Gs, Tau1Gs, Tau2Gs);
  
  Budget = KWealth + NInc - SimulateResult.tax(:,j);
  SimulateResult.c(:,j) = (Budget - SimulateResult.k(:,j+1)) / (1+TauC);

  SimulateResult.ninc(:,j) = NInc;
  SimulateResult.kinc(:,j) = KInc;
  
  SimulateResult.survival(:,j) = SurvivalRate(j);
  SimulateResult.population(:,j+1) = SimulateResult.population(:,j) .* SimulateResult.survival(:,j) / (1+PopGrowthRate);
end

% Simulate old
ProblemOld.options.Algorithm = 'brent';
ProblemOld.options.NumThreads = NumOfThreads;
ProblemOld.options.TolX = TolOpt;
%{
ValuePp = struct('form','pp','breaks',{{KGrid}},...
    'Values',[],'coefs',[],'order',[4],...
    'Method',[],'ExtrapolationOrder',[],'thread',NumOfThreads,...
    'orient','curvefit');
ValuePp.Values = reshape(permute(ValueOld, [2 1]), Jr+1, KPts);
ValuePp = myppual(ValuePp);
%}
ValuePp = tensor_pchip(KGrid, reshape(permute(ValueOld, [2 1]), Jr+1, KPts));
ProblemOld.pp = myppual(ValuePp);
ProblemOld.lb = KPrimeMin * ones(1, NumOfAgents);
ProblemOld.shared_data = [
    TauC;
    Sigma1;
    Chi;
    ];

SimulateResult.s(:,Jw+1:J) = 0;
SimulateResult.n(:,Jw+1:J) = 0;

for j=Jw+1:1:J
  % compute value interpolation vector function index
  VFutureIdx = (j-Jw+1)*ones(1, NumOfAgents);
  % use this period distribution as data
  KInc = r*(SimulateResult.k(:,j) + Tr);
  KWealth = (1+r)*(SimulateResult.k(:,j) + Tr);
  SimulateResult.tax(:,j) = GS(KInc, Tau0Gs, Tau1Gs, Tau2Gs);
  Budget = KWealth - SimulateResult.tax(:,j) + Tss;

  % construct problem
  %{
  if j==J
      ProblemOld.lb = zeros(1,NumOfAgents);
  else
      ProblemOld.lb = KPrimeMin*ones(1,NumOfAgents);
  end
  %}
  
  ProblemOld.ub = min(Budget(:)' - CMin, KPrimeMax);
  ProblemOld.x0 = (ProblemOld.lb + ProblemOld.ub) / 2;

  ProblemOld.Data = [
      Budget(:)'
	  ];
  ProblemOld.idx = VFutureIdx;
  
  [X, V] = dpopt_utility_old(ProblemOld);
  
  SimulateResult.k(:,j+1) = X;
  SimulateResult.v(:,j) = V;

  SimulateResult.c(:,j) = (Budget - SimulateResult.k(:,j+1)) / (1+TauC);
  SimulateResult.kinc(:,j) = KInc;
  SimulateResult.ninc(:,j) = 0;
  SimulateResult.survival(:,j) = SurvivalRate(j);
  SimulateResult.population(:,j+1) = SimulateResult.population(:,j) .* SimulateResult.survival(:,j) / (1+PopGrowthRate);
end

% write J+1 capital and population to empty
SimulateResult.k(:,end) = [];
SimulateResult.population(:,end) = [];

ExitFlag = 0;
end
