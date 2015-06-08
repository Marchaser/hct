function [VfiResult, ExitFlag] = VFI_NO_RISK(Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Solve household value function iteration

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
ProblemOld.options.Algorithm = 'golden';
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

% endogenous borrowing constraints
EBC = 0;

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
    
    ProblemOld.ub = Budget(:)' - CMin;
    ProblemOld.lb = max(KPrimeMin, EBC) * ones(1, KPts);
    ProblemOld.ub = max(ProblemOld.ub, ProblemOld.lb);
    ProblemOld.x0 = (ProblemOld.lb + ProblemOld.ub) / 2;
    
    % endogenous borrowing constraints grow by adding current period
    % minimum income
    EBC = (EBC-Tss) / (1+r);
    
    % solve problem
    [X, V, OptExitFlag] = dpopt_utility_old(ProblemOld);
    
    % update policy and value
    PolicyOld(:,j) = X;
    ValueOld(:,j) = V * Beta * SurvivalRate(Jw+j);
end

% young's problem
KInc = r * (KMesh + Tr);

% NOTE(wenlan): The dpopt problem takes (s,n,k') as controls. One
% nonlinear constraint is for c; one linear constraint is for l.
ProblemYoung.options.Algorithm = 'donlp2';
ProblemYoung.options.NumThreads = NumOfThreads;
ProblemYoung.options.TolX = TolOpt;
ProblemYoung.options.TolCon = TolCon;
ProblemYoung.options.ViolationLimit = 1e-6;
ProblemYoung.lb = [
    SMin * ones(1, NumOfProblems)
    NMin * ones(1, NumOfProblems)
    KPrimeMin * ones(1, NumOfProblems)
    ];

% NOTE(wenlan): I'm avoiding all extrapolations in VFI for robustness. It's
% guranteed that the domain of grid is large enough s.t. there's no
% extrapolation in simulation
ProblemYoung.ub = [
    SMax * ones(1, NumOfProblems)
    NMax * ones(1, NumOfProblems)
    KPrimeMax * ones(1, NumOfProblems)
    ];

ProblemYoung.Aineq = repmat([1;1;0], 1, NumOfProblems);
ProblemYoung.blb = zeros(1, NumOfProblems);
ProblemYoung.bub = ones(1, NumOfProblems) - LMin;
ProblemYoung.clb = CMin * ones(1, NumOfProblems);
ProblemYoung.cub = 1e20 * ones(1, NumOfProblems);

ProblemYoung.shared_data = [
    r; w; Tr;
    Tau0Gs; Tau1Gs; Tau2Gs;
    TauSs; TauC;
    Alpha; Rho;
    Chi; Sigma1; Sigma2;
    1;
    0.9955;
    1;
    ];

ProblemYoung.Data = [...
    AMesh(:)';
    KMesh(:)';
    HMesh(:)';
    KInc(:)';
    ];

% NOTE(wenlan): vector function index to integrate over future epsilons.
% Keep in mind, a goes frst in memory
ProblemYoung.idx = repmat([1:APts], EpsilonPts, KPts * HPts);

% initiate young's value function using old's
ValueYoung = zeros([APts KPts HPts Jw+1]);
PolicyYoung = zeros([3 APts KPts HPts Jw]);
ValueYoung(:,:,:,Jw+1) = repmat(reshape(ValueOld(:,1), [1 KPts 1]), [APts 1 HPts]);

% use old's solution in capital as initial guess
X = [
    1e-3 * ones(1, NumOfProblems)
    1e-3 * ones(1, NumOfProblems)
    reshape(repmat(reshape(PolicyOld(:,1), [1 KPts 1]), [APts 1 HPts]), 1, NumOfProblems)
    ];

for j = Jw:-1:1
    V0 = reshape(ValueYoung(:,:,:,j+1), [APts KPts HPts]);
    
    V0Pp = tensor_pchip({KGrid, HGrid}, V0);
    ProblemYoung.pp = myppual(V0Pp);
    
    % account for endogenous borrowing constraints
    ProblemYoung.lb(3,:) = max(KPrimeMin, EBC) * ones(1, NumOfProblems);
    ProblemYoung.ub(3,:) = max(ProblemYoung.ub(3,:),ProblemYoung.lb(3,:));
    EBC = (EBC-Tr) / (1+r);
    
    % Note(wenlan): use last iteration solution as initial guess for this
    % iteration
    ProblemYoung.x0 = X;
    switch ControlSpec
        case 'Exn'
            ProblemYoung.x0(2,:) = reshape(BenchResult.VfiResult.PolicyYoung(2,:,:,:,j), 1, []);
            ProblemYoung.lb(2,:) = ProblemYoung.x0(2,:);
            ProblemYoung.ub(2,:) = ProblemYoung.x0(2,:);
    end
    
    [X, V, OptExitFlag] = dpopt_utility(ProblemYoung);
    
    % update policy and value
    PolicyYoung(:,:,:,:,j) = reshape(X, [3 APts KPts HPts]);
    ValueYoung(:,:,:,j) = reshape(V, ProblemSize) * Beta * SurvivalRate(j);
end

% return
VfiResult = v2struct(PolicyOld, ValueOld, PolicyYoung, ValueYoung);

ExitFlag = 0;
end

function tempplot
mesh(squeeze(KMesh(1,:,:)), squeeze(HMesh(1,:,:)), squeeze(ValueOld(:,:,j)));
xlabel('k');
ylabel('h');

figure;
line(KGrid,ValueOld(:,j));
xlabel('k');
figure;
line(KGrid,PolicyOld(:,j));
xlabel('k');

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
KShrinkFactor = 2;
HShrinkFactor = 2;
KCrudeIdx = 1:KShrinkFactor:length(KGrid);
HCrudeIdx = 1:HShrinkFactor:length(HGrid);
KCrudeGrid = KGrid(KCrudeIdx);
HCrudeGrid = HGrid(HCrudeIdx);
[KCrudeMesh,HCrudeMesh] = ndgrid(KCrudeGrid,HCrudeGrid);
PolicyYoungCrude = reshape(PolicyYoung(2,1,KCrudeIdx,HCrudeIdx,j),1,length(KCrudeIdx),length(HCrudeIdx));
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

mesh(squeeze(KMesh(1,:,:)), squeeze(HMesh(1,:,:)), squeeze(PolicyYoung(2,1,:,:,j)));
xlabel('k');
ylabel('h');

mesh(squeeze(KMesh(1,:,:)), squeeze(HMesh(1,:,:)), squeeze(PolicyYoung(3,1,:,:,j)));
xlabel('k');
ylabel('h');
end
