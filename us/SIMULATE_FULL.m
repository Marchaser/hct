function [SimulateResult, ExitFlag] = SIMULATE(VfiResult, Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Simulate economy, return distribution

% unpack structure
v2struct(Params);
v2struct(VfiResult);

% write TauGs to empty on factor tax
switch TaxSpec
    case 'Factor'
        Tau0Gs = [];
        Tau1Gs = [];
        Tau2Gs = [];
end

% initialize storage
% Note(wenlan): variable starting with lower case is individual variable,
% to be distinguished from the aggregate variable starting with upper case
SimulateResult.c = zeros(NumOfAgents, J);
SimulateResult.k = zeros(NumOfAgents, J+1);
SimulateResult.h = zeros(NumOfAgents, J);
SimulateResult.s = zeros(NumOfAgents, J);
SimulateResult.n = zeros(NumOfAgents, J);
SimulateResult.l = zeros(NumOfAgents, J);

SimulateResult.v = zeros(NumOfAgents, J);

SimulateResult.kinc = zeros(NumOfAgents, J);
SimulateResult.ninc = zeros(NumOfAgents, J);
SimulateResult.ntax = zeros(NumOfAgents, J);
SimulateResult.ktax = zeros(NumOfAgents, J);
SimulateResult.tax = zeros(NumOfAgents, J);
SimulateResult.tss = zeros(NumOfAgents, J);
SimulateResult.survival = ones(NumOfAgents, J);
SimulateResult.population = ones(NumOfAgents, J+1);

% initial distribution
SimulateResult.k(:,1) = K1;
SimulateResult.h(:,1) = H1;

% simulate young
% NOTE(wenlan): Construct policy function spline. Keep in mind,
% PolicyYoung has the memory order (control,A,K,H,Age), permute K and H to
% last two dimensions for construction
switch(ControlSpec)
    case {'Base','Exn'}
        PolicyPp = tensor_pchip({KGrid,HGrid}, reshape(permute(PolicyYoung, [1 5 2 3 4]), [], KPts, HPts));
        PolicyVecSize = [3 Jw APts];
        ValuePp = tensor_pchip({KGrid,HGrid}, ValueYoung(:,:,:,1));
    case 'Exhc'
        % (control,A,K,KBar,HBar,Age)
        PolicyPp = tensor_pchip({KGrid,KBarGrid,HBarGrid}, reshape(permute(PolicyYoung, [1 6 2 3 4 5]), [], KPts, KBarPts, HBarPts));
        PolicyVecSize = [2 Jw APts];
        ValuePp = tensor_pchip({KGrid,KBarGrid,HBarGrid}, ValueYoung(:,:,:,:,1));
end
PolicyPp = myppual(PolicyPp);
ValuePp = myppual(ValuePp);

%{
switch(ControlSpec)
    case {'Base','Exn'}
        SimulateResult.v = myppual(ValuePp, [SimulateResult.k(:,1)'; SimulateResult.h(:,1)'], [], AIdx(:)');
    case 'Exhc'
        SimulateResult.v = myppual(ValuePp, [SimulateResult.k(:,1)'; BenchSimulationK(:,1)'; BenchSimulationH(:,1)'], [], AIdx(:)');
end
%}


for j=1:1:Jw
    switch(ControlSpec)
        case {'Base','Exn'}
            ValuePp = tensor_pchip({KGrid,HGrid}, ValueYoung(:,:,:,j));
            ValuePp = myppual(ValuePp);
            SimulateResult.v(:,j) = myppual(ValuePp, [SimulateResult.k(:,j)'; SimulateResult.h(:,j)'], [], AIdx(:)');
        case 'Exhc'
            SimulateResult.v = myppual(ValuePp, [SimulateResult.k(:,1)'; BenchSimulationK(:,1)'; BenchSimulationH(:,1)'], [], AIdx(:)');
    end
    
    % interpolate to get policy
    switch(ControlSpec)
        case 'Base'
            PolicyInterpIdx = sub2ind(PolicyVecSize, ...
                repmat([1:3]', 1, NumOfAgents), ...
                j*ones(3, NumOfAgents), ...
                repmat(AIdx(:)', 3, 1));
            X = myppual(PolicyPp, [SimulateResult.k(:,j)'; SimulateResult.h(:,j)'], [], PolicyInterpIdx);
            % extract policy
            SimulateResult.s(:,j) = X(1,:);
            SimulateResult.n(:,j) = X(2,:);
            SimulateResult.k(:,j+1) = X(3,:);
        case 'Exn'
            PolicyInterpIdx = sub2ind(PolicyVecSize, ...
                repmat([1:3]', 1, NumOfAgents), ...
                j*ones(3, NumOfAgents), ...
                repmat(AIdx(:)', 3, 1));
            X = myppual(PolicyPp, [SimulateResult.k(:,j)'; SimulateResult.h(:,j)'], [], PolicyInterpIdx);
            % extract policy
            SimulateResult.s(:,j) = X(1,:);
            SimulateResult.n(:,j) = BenchSimulationN(:,j);
            SimulateResult.k(:,j+1) = X(3,:);
        case 'Exhc'
            PolicyInterpIdx = sub2ind(PolicyVecSize, ...
                repmat([1:2]', 1, NumOfAgents), ...
                j*ones(2, NumOfAgents), ...
                repmat(AIdx(:)', 2, 1));
            X = myppual(PolicyPp, [SimulateResult.k(:,j)';BenchSimulationK(:,j)';BenchSimulationH(:,j)'], [], PolicyInterpIdx);
            % extract policy
            SimulateResult.s(:,j) = BenchSimulationS(:,j);
            SimulateResult.n(:,j) = X(1,:);
            SimulateResult.k(:,j+1) = X(2,:);
    end
    
    % enforce boundary
    SimulateResult.s(:,j) = max(SimulateResult.s(:,j), 0);
    SimulateResult.n(:,j) = max(SimulateResult.n(:,j), 0);
    
    % next period h
    SimulateResult.h(:,j+1) = EpsilonShock(:,j) .* ...
        ((1-Rho)*SimulateResult.h(:,j) + AShock(:) .* (SimulateResult.s(:,j).*SimulateResult.h(:,j)).^Alpha);
    
    NInc = w * (1-TauSs) * SimulateResult.n(:,j) .* SimulateResult.h(:,j);
    SimulateResult.tss(:,j) = w*TauSs*SimulateResult.n(:,j) .* SimulateResult.h(:,j);
    KInc = r * (SimulateResult.k(:,j) + Tr);
    KWealth = (1+r) * (SimulateResult.k(:,j) + Tr);
    
    switch TaxSpec
        case 'Total'
            SimulateResult.tax(:,j) = GS(KInc + NInc, Tau0Gs, Tau1Gs, Tau2Gs);
        case 'Factor'
            SimulateResult.ktax(:,j) = GS(KInc, Tau0K, Tau1K, Tau2K);
            SimulateResult.ntax(:,j) = GS(NInc, Tau0N, Tau1N, Tau2N);
            SimulateResult.tax(:,j) = SimulateResult.ktax(:,j) + SimulateResult.ntax(:,j);
    end
    
    Budget = KWealth + NInc - SimulateResult.tax(:,j);
    SimulateResult.c(:,j) = (Budget - SimulateResult.k(:,j+1)) / (1+TauC);
    SimulateResult.l(:,j) = 1-SimulateResult.s(:,j)-SimulateResult.n(:,j);
    
    % enforce boundary
    SimulateResult.c(:,j) = max(SimulateResult.c(:,j), 0);
    SimulateResult.l(:,j) = max(SimulateResult.l(:,j), 0);
    
    SimulateResult.ninc(:,j) = NInc;
    SimulateResult.kinc(:,j) = KInc;
    
    SimulateResult.survival(:,j) = SurvivalRate(j);
    SimulateResult.population(:,j+1) = SimulateResult.population(:,j) .* SimulateResult.survival(:,j) / (1+PopGrowthRate);
end

% Simulate old
PolicyPp = tensor_pchip(KGrid, reshape(permute(PolicyOld, [2 1]), Jr, KPts));
SimulateResult.s(:,Jw+1:J) = 0;
SimulateResult.n(:,Jw+1:J) = 0;
SimulateResult.l(:,Jw+1:J) = 1;

for j=Jw+1:1:J
    ValuePp = tensor_pchip(KGrid, ValueOld(:,j-Jw));
    ValuePp = myppual(ValuePp);
    SimulateResult.v(:,j) = myppual(ValuePp, [SimulateResult.k(:,j)'], [], ones(1,NumOfAgents));
            
    % compute value interpolation vector function index
    PolicyInterpIdx = (j-Jw)*ones(1, NumOfAgents);
    % use this period distribution as data
    KInc = r*(SimulateResult.k(:,j) + Tr);
    KWealth = (1+r)*(SimulateResult.k(:,j) + Tr);
    
    switch TaxSpec
        case 'Total'
            SimulateResult.tax(:,j) = GS(KInc, Tau0Gs, Tau1Gs, Tau2Gs);
        case 'Factor'
            SimulateResult.ktax(:,j) = GS(KInc, Tau0K, Tau1K, Tau2K);
            SimulateResult.ntax(:,j) = 0;
            SimulateResult.tax(:,j) = SimulateResult.ktax(:,j) + SimulateResult.ntax(:,j);
    end
    
    Budget = KWealth - SimulateResult.tax(:,j) + Tss;
    
    X = myppual(PolicyPp, [SimulateResult.k(:,j)'], [], PolicyInterpIdx);
    
    SimulateResult.k(:,j+1) = X;
    
    SimulateResult.c(:,j) = (Budget - SimulateResult.k(:,j+1)) / (1+TauC);
    
    % enforce boundary
    SimulateResult.c(:,j) = max(SimulateResult.c(:,j), 0);
    
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
