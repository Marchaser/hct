function [VfiResult, ExitFlag] = VFI_EXHC(Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Solve household value function iteration
% with exogenous human capital

% unpack structures
v2struct(Params);

if (r < 0)
    VfiResult = [];
    ExitFlag = -1;
    return;
end

% old's problem
Budget = (1+r)*(KGrid+Tr) - GS(r*(KGrid+Tr), Tau0Gs, Tau1Gs, Tau2Gs) + Tss;

% construct problem
ProblemOld.options.Algorithm = 'brent';
ProblemOld.options.NumThreads = 1;
ProblemOld.options.TolX = TolOpt;
ProblemOld.lb = KPrimeMin * ones(1, KPts);
ProblemOld.ub = Budget(:)' - CMin;
ProblemOld.x0 = (ProblemOld.lb + ProblemOld.ub) / 2;

ProblemOld.shared_data = [
    TauC
    Sigma1
    Chi
    ];

ProblemOld.Data = [
	Budget(:)'
	];

ProblemOld.idx = ones(1, KPts);

% initiate value function
ValueOld = zeros([KPts Jr+1]);
PolicyOld = zeros([KPts Jr]);

% solve old's problem
for j = Jr:-1:1
    % generate v0 from next period
    V0 = ValueOld(:, j+1);
    %{
    V0Pp = struct('form','pp','breaks',{{KGrid}},...
        'Values',V0,'coefs',[],'order',[4],...
        'Method',[],'ExtrapolationOrder',[],'thread',1,...
        'orient','curvefit');
    V0Pp = myppual(V0Pp);
    %}
    V0Pp = tensor_pchip(KGrid, V0);
    ProblemOld.pp = myppual(V0Pp);
    
    if (j==Jr)
        ProblemOld.ub = zeros(1, KPts);
        ProblemOld.lb = zeros(1, KPts);
        ProblemOld.x0 = (ProblemOld.lb + ProblemOld.ub) / 2;
    else
        ProblemOld.ub = Budget(:)' - CMin;
        ProblemOld.lb = KPrimeMin * ones(1, KPts);
        ProblemOld.x0 = (ProblemOld.lb + ProblemOld.ub) / 2;
    end
    
    % solve problem
    [X, V] = dpopt_utility_old(ProblemOld);
    
    % update policy and value
    PolicyOld(:,j) = X;
    ValueOld(:,j) = V * Beta * SurvivalRate(Jw+j);
end

% young's problem
[AMesh, KBarMesh, HBarMesh] = ndgrid(AGrid, KBarGrid, HBarGrid);
NumOfProblems = APts * KPts * KBarPts * HBarPts;
KInc = r * (KGrid + Tr);

% NOTE(wenlan): The dpopt problem takes (n,k') as controls. One
% nonlinear constraint is for c.
ProblemYoung.options.Algorithm = 'donlp2';
ProblemYoung.options.NumThreads = NumOfThreads;
ProblemYoung.options.TolX = TolOpt;
ProblemYoung.options.TolCon = TolCon;
ProblemYoung.options.ViolationLimit = 1e-6;
ProblemYoung.lb = [
    NMin * ones(1, NumOfProblems)
    KPrimeMin * ones(1, NumOfProblems)
    ];

% NOTE(wenlan): I'm avoiding all extrapolations in VFI for robustness. It's
% guranteed that the domain of grid is large enough s.t. there's no
% extrapolation in simulation
ProblemYoung.ub = [
    NMax * ones(1, NumOfProblems)
    KPrimeMax * ones(1, NumOfProblems)
    ];

ProblemYoung.clb = CMin * ones(1, NumOfProblems);
ProblemYoung.cub = 1e20 * ones(1, NumOfProblems);

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

% NOTE(wenlan): vector function index to integrate over future epsilons.
% Keep in mind, a goes frst in memory
ProblemYoung.idx = repmat([1:APts], EpsilonPts, KPts * KBarPts * HBarPts);

% initiate young's value function using old's
ValueYoung = zeros([APts KPts KBarPts HBarPts Jw+1]);
PolicyYoung = zeros([2 APts KPts KBarPts HBarPts Jw]);
ValueYoung(:,:,:,:,Jw+1) = repmat(reshape(ValueOld(:,1), [1 KPts 1 1]), [APts 1 KBarPts HBarPts]);

% use old's solution in capital as initial guess
X = [
    1e-3 * ones(1, APts*KPts*KBarPts*HBarPts)
    reshape(repmat(reshape(PolicyOld(:,1), [1 KPts 1 1]), [APts 1 KBarPts HBarPts]), 1, APts*KPts*KBarPts*HBarPts)
    ];
X = reshape(X, [2 APts KPts KBarPts HBarPts]);

for j = Jw:-1:1
  V0 = reshape(ValueYoung(:,:,:,:,j+1), [APts KPts KBarPts HBarPts]);
  
  % construct spline
  %{
  V0Pp = struct('form','pp','breaks',{{KGrid,HGrid}},...
    'Values',V0,'coefs',[],'order',[4 4],...
    'Method',[],'ExtrapolationOrder',[],'thread',NumOfThreads,...
    'orient','curvefit');
  V0Pp = myppual(V0Pp);
  %}
  
  %{
    NOTE(wenlan): Here's what I'm doing.
  (1) V0 has the memory order (A,K,KBar,HBar). The reason for organizing V0
  in this order is I can interpolate over the last two psuedo states.
  (2) I treat (A,K) as vector function and interpolate over (KBarPrime,HBarPrime).
  (3) Then I permute the interpolated results to order (A,KBarPrime,HBarPrime,K).
  (4) Next I treat (A,KBarPrime,HBarPrime) as vector function and construct
  splines over K.
  (5) I pass the interpolation over K to optimizers
  %}
  
  % Do a construction w.r.t. KBar and HBar, treating (A K) as vector functions
  V0PpAtKBarHBarPp = tensor_pchip({KBarGrid, HBarGrid}, reshape(V0, [], KBarPts, HBarPts));
  V0PpAtKBarHBarPp = myppual(V0PpAtKBarHBarPp);
  
  KBarPrime = BenchKpPolicy(j,:,:,:);
  KBarPrime = repmat(KBarPrime(:)', EpsilonPts*KPts, 1);
  % Change memory order to Epsilon A K KBar HBar
  KBarPrime = permute(reshape(KBarPrime, [EpsilonPts KPts APts KBarPts HBarPts]), [1 3 2 4 5]);
  
  SBar = squeeze(BenchSPolicy(j,:,:,:));
  HBarPrimeBeforeShock = (1-Rho)*HBarMesh + AMesh.*(SBar.*HBarMesh).^Alpha;
  HBarPrimeBeforeShock = repmat(HBarPrimeBeforeShock(:)', EpsilonPts*KPts, 1);
  % Change memory order to Epsilon A K KBar HBar
  HBarPrimeBeforeShock = permute(reshape(HBarPrimeBeforeShock, [EpsilonPts, KPts APts KBarPts HBarPts]), [1 3 2 4 5]);
  % Expand to accomodate Epsilon shock
  HBarPrime = HBarPrimeBeforeShock .* repmat(EpsilonGrid(:), [1 APts KPts KBarPts HBarPts]);
  
  % Interpolate over corresponding A K idx
  V0InterpIdx = repmat(reshape(1:APts*KPts, [1 APts KPts 1 1]), [EpsilonPts 1 1 KBarPts HBarPts]);
  
  % interpolate over all combinitions of KBarPrime and HBarPrime
  V0AtKBarPrimeHBarPrime = myppual(V0PpAtKBarHBarPp, [KBarPrime(:)';HBarPrime(:)'], [], V0InterpIdx(:)');
  
  % Integaret over Epsilon shock
  EV0AtKBarPrimeHBarPrime = EpsilonP * reshape(V0AtKBarPrimeHBarPrime, EpsilonPts, []);
  
  % Change memroy order to A KBar HBar K
  EV0AtKBarPrimeHBarPrime = reshape(EV0AtKBarPrimeHBarPrime, [APts KPts KBarPts HBarPts]);
  EV0AtKBarPrimeHBarPrime = permute(EV0AtKBarPrimeHBarPrime, [1 3 4 2]);
  
  % Construct pp w.r.t. K
  V0KPp = tensor_pchip(KGrid, reshape(EV0AtKBarPrimeHBarPrime, [], KPts));
  
  % ValueYoung's memory order is A K KBar HBar
  % I first construct index in order K A KBar HBar, then permute
  V0KInterpIdx = repmat(1:APts*KBarPts*HBarPts, KPts, 1);
  V0KInterpIdx = permute(reshape(V0KInterpIdx, [KPts APts KBarPts HBarPts]), [2 1 3 4]);
  ProblemYoung.idx = V0KInterpIdx(:)';
  
  % Similarly, deal with SBarData and HBarData
  SBarData = repmat(reshape(BenchSPolicy(j,:,:,:),1,[]), KPts, 1);
  SBarData = permute(reshape(SBarData, [KPts APts KBarPts, HBarPts]), [2 1 3 4]);
  HBarData = repmat(reshape(HBarMesh(:,:,:),1,[]), KPts, 1);
  HBarData = permute(reshape(HBarData, [KPts APts KBarPts, HBarPts]), [2 1 3 4]);
  
  % deal with KGrid and KInc
  KData = repmat(reshape(KGrid, [KPts 1 1 1]), [1 APts KBarPts HBarPts]);
  KData = permute(KData, [2 1 3 4]);
  KIncData = repmat(reshape(KInc, [KPts 1 1 1]), [1 APts KBarPts HBarPts]);
  KIncData = permute(KIncData, [2 1 3 4]);
  
  % Pass to dpopt
  ProblemYoung.pp = myppual(V0KPp);
  ProblemYoung.x0 = X;
  ProblemYoung.Data = [...
      KData(:)';
      SBarData(:)';
      HBarData(:)';
      KIncData(:)';
      ];
  [X, V, OptExitFlag] = dpopt_utility_exhc(ProblemYoung);
 
  % update policy and value
  PolicyYoung(:,:,:,:,:,j) = reshape(X, [2 APts KPts KBarPts HBarPts]);
  ValueYoung(:,:,:,:,j) = reshape(V, [APts KPts KBarPts HBarPts]) * Beta * SurvivalRate(j);
end

% return
VfiResult = v2struct(PolicyOld, ValueOld, PolicyYoung, ValueYoung);

ExitFlag = 0;
end

function tempplot
mesh(squeeze(KMesh(1,:,:)), squeeze(HMesh(1,:,:)), squeeze(ValueOld(:,:,j)));
xlabel('k');
ylabel('h');

mesh(squeeze(KMesh(1,:,:)), squeeze(HMesh(1,:,:)), squeeze(ValueYoung(1,:,:,j+1)));
xlabel('k');
ylabel('h');
% NOTE(wenlan): Interpolate over fine grids to observe behavior of
% interpolation
KFineGrid = linspace(KMin,KMax,100);
HFineGrid = linspace(HMin,HMax,100);
[KFineMesh,HFineMesh] = ndgrid(KFineGrid,HFineGrid);
% Supply as vectorized data to myppual for interpolation, interpolate over
% the first vector function
VFutureFine = myppual(ProblemYoung.pp,[KFineMesh(:)';HFineMesh(:)'],[],ones(1,length(KFineMesh(:))));
% Plot the interpolated future value over fine meshes
mesh(KFineMesh, HFineMesh, reshape(VFutureFine,100,100));

mesh(squeeze(KMesh(1,:,:)), squeeze(HMesh(1,:,:)), squeeze(ValueYoung(1,:,:,j)));
xlabel('k');
ylabel('h');

mesh(squeeze(KMesh(:,:,10)), squeeze(AMesh(:,:,10)), squeeze(ValueYoung(:,:,10,j)));
xlabel('k');
ylabel('a');

figure;
mesh(squeeze(KMesh(1,:,:)), squeeze(HMesh(1,:,:)), squeeze(PolicyYoung(1,1,:,:,j)));
xlabel('k');
ylabel('h');
% Note(wenlan): This is to test how many grids to produce satisfactory
% interpolation of policy function
KShrinkFactor = 3;
HShrinkFactor = 3;
KCrudeIdx = 1:KShrinkFactor:length(KGrid);
HCrudeIdx = 1:HShrinkFactor:length(HGrid);
KCrudeGrid = KGrid(KCrudeIdx);
HCrudeGrid = HGrid(HCrudeIdx);
[KCrudeMesh,HCrudeMesh] = ndgrid(KCrudeGrid,HCrudeGrid);
PolicyYoungCrude = reshape(PolicyYoung(1,1,KCrudeIdx,HCrudeIdx,j),1,length(KCrudeIdx),length(HCrudeIdx));
figure;
mesh(KCrudeMesh,HCrudeMesh,squeeze(PolicyYoungCrude))
PolicyYoungCrudePp = tensor_pchip({KCrudeGrid,HCrudeGrid},PolicyYoungCrude);
KFineGrid = linspace(KMin,KMax,100);
HFineGrid = linspace(HMin,HMax,100);
[KFineMesh,HFineMesh] = ndgrid(KFineGrid,HFineGrid);
PolicyYoungInterpFine = myppual(PolicyYoungCrudePp,[KFineMesh(:)';HFineMesh(:)'],[],ones(1,length(KFineMesh(:))));
figure;
mesh(KFineMesh, HFineMesh, reshape(PolicyYoungInterpFine,100,100));

mesh(squeeze(KMesh(:,:,1)), squeeze(AMesh(:,:,1)), squeeze(PolicyYoung(1,:,:,1,j)));
xlabel('k');
ylabel('a');

figure;
mesh(squeeze(KBarMesh(1,:,:)), squeeze(HBarMesh(1,:,:)), squeeze(ValueYoung(1,20,:,:,j)));
xlabel('k');
ylabel('h');

mesh(squeeze(KBarMesh(1,:,:)), squeeze(HBarMesh(1,:,:)), squeeze(PolicyYoung(2,1,1,:,:,j)));
xlabel('k');
ylabel('h');

mesh(squeeze(KMesh(1,:,:)), squeeze(HMesh(1,:,:)), squeeze(PolicyYoung(3,1,:,:,j)));
xlabel('k');
ylabel('h');
end
