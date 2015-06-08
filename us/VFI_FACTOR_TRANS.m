function [VfiResult, ExitFlag] = VFI_FACTOR_TRANS(Params, VfiSsResult, Rt, Wt, Trt, Tsst)
% Taxing Human Capital
% @Author: Wenlan Luo
% Solve household value function iteration, with factor income tax code

r = [];
w = [];

% unpack structures
v2struct(Params);
% num of transition periods
% TTrans = length(Rt);
% extend to cover steady state
% Rt = [Rt Rt(end)*ones(1,J)];
% Wt = [Wt Wt(end)*ones(1,J)];
% Trt = [Trt Trt(end)*ones(1,J)];
% Tsst = [Tsst Tsst(end)*ones(1,J)];

% TauGs should never be used. Write it to empty to avoid mistake
Tau0Gs = [];
Tau1Gs = [];
Tau2Gs = [];

% construct problem
ProblemOld.options.Algorithm = 'brent';
ProblemOld.options.NumThreads = 1;
ProblemOld.options.TolX = TolOpt;

ProblemOld.shared_data = [
    TauC
    Sigma1
    Chi
    ];

ProblemOld.idx = ones(1, KPts);

% initiate value function
ValueOld = zeros([KPts Jr+1 TTrans+1]);
ValueOld(:,:,TTrans+1) = VfiSsResult.ValueOld;
PolicyOld = zeros([KPts Jr TTrans]);

% each periods
for t = TTrans:-1:1
    % solve old's problem
    for j = Jr:-1:1
        % generate v0 from next period
        V0 = ValueOld(:, j+1, t+1);
        
        V0Pp = tensor_pchip(KGrid, V0);
        ProblemOld.pp = myppual(V0Pp);
        
        % old's problem
        Budget = (1+Rt(t))*(KGrid+Trt(t)) - GS(Rt(t)*(KGrid+Trt(t)), Tau0K, Tau1K, Tau2K) + Tsst(t);
        ProblemOld.Data = [
            Budget(:)'
            ];
        ProblemOld.lb = KPrimeMin * ones(1, KPts);
        ProblemOld.ub = Budget(:)' - CMin;
        ProblemOld.x0 = (ProblemOld.lb + ProblemOld.ub) / 2;
        
        % solve problem
        [X, V] = dpopt_utility_old(ProblemOld);
        
        % update policy and value
        PolicyOld(:,j,t) = X;
        ValueOld(:,j,t) = V * Beta * SurvivalRate(Jw+j);
    end
end

% young's problem

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

% NOTE(wenlan): I'm avoiding all extrapolations in VFI for robustness. It's
% guranteed that the domain of grid is large enough s.t. there's no
% extrapolation in simulation


% NOTE(wenlan): vector function index to integrate over future epsilons.
% Keep in mind, a goes frst in memory
ProblemYoung.idx = repmat([1:APts], EpsilonPts, KPts * HPts);

% initiate young's value function using old's
ValueYoung = zeros([APts KPts HPts Jw+1 TTrans+1]);
ValueYoung(:,:,:,:,TTrans+1) = VfiSsResult.ValueYoung;
PolicyYoung = zeros([3 APts KPts HPts Jw TTrans]);
for t=TTrans:-1:1
    %
    ValueYoung(:,:,:,Jw+1,t+1) = repmat(reshape(ValueOld(:,1,t+1), [1 KPts 1]), [APts 1 HPts]);
    % use steady state solution as initial guess
    X = reshape(VfiSsResult.PolicyYoung(:,:,:,:,Jw), 3, []);
    
    for j = Jw:-1:1
        V0 = reshape(ValueYoung(:,:,:,j+1,t+1), [APts KPts HPts]);
        
        % construct spline
        V0Pp = tensor_pchip({KGrid, HGrid}, V0);
        ProblemYoung.pp = myppual(V0Pp);
        
        % Note(wenlan): use last iteration solution as initial guess for this
        % iteration
        ProblemYoung.x0 = X;
        
        ProblemYoung.shared_data = [
            Rt(t);Wt(t);Trt(t);
            Tau0N; Tau1N; Tau2N;
            TauSs; TauC;
            Alpha; Rho;
            Chi; Sigma1; Sigma2;
            EpsilonPts;
            EpsilonGrid(:);
            EpsilonP(:);
            ];
        
        KInc = Rt(t) * (KMesh + Trt(t));
        KTax = GS(KInc, Tau0K, Tau1K, Tau2K);
        ProblemYoung.Data = [...
            AMesh(:)';
            KMesh(:)';
            HMesh(:)';
            KInc(:)'-KTax(:)';
            ];
        [X, V, OptExitFlag] = dpopt_utility_factor(ProblemYoung);
        
        % update policy and value
        PolicyYoung(:,:,:,:,j,t) = reshape(X, [3 APts KPts HPts]);
        ValueYoung(:,:,:,j,t) = reshape(V, ProblemSize) * Beta * SurvivalRate(j);
    end
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
KGridFine = linspace(KMin,KMax,1000);
HGridFine = linspace(HMin,HMax,1000);
[KMeshFine,HMeshFine] = ndgrid(KGridFine,HGridFine);
% Supply as vectorized data to myppual for interpolation, interpolate over
% the first vector function
VFutureFine = myppual(ProblemYoung.pp,[KMeshFine(:)';HMeshFine(:)'],[],ones(1,length(KMeshFine(:))));
% Plot the interpolated future value over fine meshes
mesh(KMeshFine, HMeshFine, reshape(VFutureFine,1000,1000));

mesh(squeeze(KMesh(1,:,:)), squeeze(HMesh(1,:,:)), squeeze(ValueYoung(1,:,:,j)));
xlabel('k');
ylabel('h');

mesh(squeeze(KMesh(:,:,10)), squeeze(AMesh(:,:,10)), squeeze(ValueYoung(:,:,10,j)));
xlabel('k');
ylabel('a');

mesh(squeeze(KMesh(1,:,:)), squeeze(HMesh(1,:,:)), squeeze(PolicyYoung(1,1,:,:,j)));
xlabel('k');
ylabel('h');

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
