function [SimulateResult, ExitFlag] = SIMULATE_TRANS_FULL(VfiTransResult, SimulateSsResult, Params, Rt, Wt, Trt, Tsst)
% Taxing Human Capital
% @Author: Wenlan Luo
% Simulate economy, return distribution

% unpack structure
v2struct(Params);
v2struct(VfiTransResult);
r = [];
w = [];
Tr = [];
Tss = [];

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

% start with steady state distribution
SimulateResult = SimulateSsResult;

% simulate young
% NOTE(wenlan): Construct policy function spline. Keep in mind,
% PolicyYoung has the memory order (control,A,K,H,Age,Cohort), permute K and H to
% last two dimensions for construction
%{
switch(ControlSpec)
    case {'Base','Exn'}
        PolicyPp = tensor_pchip({KGrid,HGrid}, reshape(permute(PolicyYoung, [1 5 6 2 3 4]), [], KPts, HPts));
        PolicyVecSize = [3 Jw TTrans APts];
        ValuePp = tensor_pchip({KGrid,HGrid}, reshape(permute(ValueYoung, [4 5 1 2 3]), [], KPts, HPts));
    case 'Exhc'
        % (control,A,K,KBar,HBar,Age)
        PolicyPp = tensor_pchip({KGrid,KBarGrid,HBarGrid}, reshape(permute(PolicyYoung, [1 6 2 3 4 5]), [], KPts, KBarPts, HBarPts));
        PolicyVecSize = [2 Jw APts];
        ValuePp = tensor_pchip({KGrid,KBarGrid,HBarGrid}, ValueYoung(:,:,:,:,1));
end
PolicyPp = myppual(PolicyPp);
ValuePp = myppual(ValuePp);

PolicyOldPp = tensor_pchip(KGrid, reshape(permute(PolicyOld, [2 3 1]), [], KPts));
PolicyOldPp = myppual(PolicyOldPp);
PolicyOldVecSize = [Jr TTrans];
%}

%{
switch(ControlSpec)
    case {'Base','Exn'}
        SimulateResult.v = myppual(ValuePp, [SimulateResult.k(:,1)'; SimulateResult.h(:,1)'], [], AIdx(:)');
    case 'Exhc'
        SimulateResult.v = myppual(ValuePp, [SimulateResult.k(:,1)'; BenchSimulationK(:,1)'; BenchSimulationH(:,1)'], [], AIdx(:)');
end
%}

% need to track state along transition path
vFull = zeros(NumOfAgents,J,TTrans);
kFull = zeros(NumOfAgents,J,TTrans);
hFull = zeros(NumOfAgents,J,TTrans);

% Simulate old
SimulateResult.s(:,Jw+1:J) = 0;
SimulateResult.n(:,Jw+1:J) = 0;
SimulateResult.l(:,Jw+1:J) = 1;

for t=1:TTrans
    switch(ControlSpec)
        case {'Base','Exn'}
            PolicyPp = tensor_pchip({KGrid,HGrid}, reshape(permute(squeeze(PolicyYoung(:,:,:,:,:,t)), [1 5 2 3 4]), [], KPts, HPts));
            PolicyVecSize = [3 Jw APts];
            ValuePp = tensor_pchip({KGrid,HGrid}, reshape(permute(squeeze(ValueYoung(:,:,:,:,t)), [4 1 2 3]), [], KPts, HPts));
            ValueVecSize = [Jw+1 APts];
        case 'Exhc'
            % (control,A,K,KBar,HBar,Age)
            PolicyPp = tensor_pchip({KGrid,KBarGrid,HBarGrid}, reshape(permute(PolicyYoung(:,:,:,:,:,:,t), [1 6 2 3 4 5]), [], KPts, KBarPts, HBarPts));
            PolicyVecSize = [2 Jw APts];
            ValuePp = tensor_pchip({KGrid,KBarGrid,HBarGrid}, reshape(permute(squeeze(ValueYoung(:,:,:,:,:,t)), [5 1 2 3 4]), [], KPts, KBarPts, HBarPts));
            ValueVecSize = [Jw+1 APts];
    end
    PolicyPp = myppual(PolicyPp);
    PolicyPp.thread = NumOfThreads;
    
    ValuePp = myppual(ValuePp);
    ValuePp.thread = NumOfThreads;
    
    
    PolicyOldPp = tensor_pchip(KGrid, reshape(permute(squeeze(PolicyOld(:,:,t)), [2 1]), [], KPts));
    PolicyOldPp = myppual(PolicyOldPp);
    PolicyOldPp.thread = NumOfThreads;
    
    ValueOldPp = tensor_pchip(KGrid, reshape(permute(squeeze(ValueOld(:,:,t)), [2 1]), [], KPts));
    ValueOldPp = myppual(ValueOldPp);
    ValueOldPp.thread = NumOfThreads;
    
    for j=J:-1:Jw+1
        % compute value interpolation vector function index
%         PolicyInterpIdx = sub2ind(PolicyOldVecSize, (j-Jw)*ones(1, NumOfAgents),t*ones(1,NumOfAgents));
        PolicyInterpIdx = (j-Jw)*ones(1,NumOfAgents);
        % use this period distribution as data
        KInc = Rt(t)*(SimulateResult.k(:,j) + Trt(t));
        KWealth = (1+Rt(t))*(SimulateResult.k(:,j) + Trt(t));
        
        SimulateResult.ktax(:,j) = GS(KInc, Tau0K, Tau1K, Tau2K);
        SimulateResult.ntax(:,j) = 0;
        SimulateResult.tax(:,j) = SimulateResult.ktax(:,j) + SimulateResult.ntax(:,j);
        
        Budget = KWealth - SimulateResult.tax(:,j) + Tsst(t);
        
        % record value
        vFull(:,j,t) = myppual(ValueOldPp, [SimulateResult.k(:,j)'], [], PolicyInterpIdx);
        kFull(:,j,t) = SimulateResult.k(:,j);
        hFull(:,j,t) = SimulateResult.h(:,j);
        
        X = myppual(PolicyOldPp, [SimulateResult.k(:,j)'], [], PolicyInterpIdx);
        
        SimulateResult.k(:,j+1) = X;
        
        SimulateResult.c(:,j) = (Budget - SimulateResult.k(:,j+1)) / (1+TauC);
        
        % enforce boundary
        SimulateResult.c(:,j) = max(SimulateResult.c(:,j), 0);
        
        SimulateResult.kinc(:,j) = KInc;
        SimulateResult.ninc(:,j) = 0;
        SimulateResult.survival(:,j) = SurvivalRate(j);
        SimulateResult.population(:,j+1) = SimulateResult.population(:,j) .* SimulateResult.survival(:,j) / (1+PopGrowthRate);
    end
    
    for j=Jw:-1:1
        %
        ValueInterpIdx = sub2ind(ValueVecSize, ...
            j*ones(1,NumOfAgents), ...
            AIdx(:)');
        
        switch(ControlSpec)
            case {'Base','Exn'}
                vFull(:,j,t) = myppual(ValuePp, [SimulateResult.k(:,j)';SimulateResult.h(:,j)'], [], ValueInterpIdx);
            case 'Exhc'
                vFull(:,j,t) = myppual(ValuePp, [SimulateResult.k(:,j)';BenchSimulationK(:,j)';BenchSimulationH(:,j)'], [], ValueInterpIdx);
        end
        
        kFull(:,j,t) = SimulateResult.k(:,j);
        hFull(:,j,t) = SimulateResult.h(:,j);
        
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
        
        NInc = Wt(t) * (1-TauSs) * SimulateResult.n(:,j) .* SimulateResult.h(:,j);
        SimulateResult.tss(:,j) = Wt(t)*TauSs*SimulateResult.n(:,j) .* SimulateResult.h(:,j);
        KInc = Rt(t) * (SimulateResult.k(:,j) + Trt(t));
        KWealth = (1+Rt(t)) * (SimulateResult.k(:,j) + Trt(t));
        
        SimulateResult.ktax(:,j) = GS(KInc, Tau0K, Tau1K, Tau2K);
        SimulateResult.ntax(:,j) = GS(NInc, Tau0N, Tau1N, Tau2N);
        SimulateResult.tax(:,j) = SimulateResult.ktax(:,j) + SimulateResult.ntax(:,j);
        
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
    % aggregation takes place
    SimulateResult.k(:,end) = [];
    SimulateResult.population(:,end) = [];
    
    v2struct(SimulateResult);
    C(t) = WMean(c(:), population(:));
    TaxC(t) = C(t) * TauC;
    K(t) = WMean(k(:), population(:));
    L(t) = WMean(h(:).*n(:), population(:));
    Tss(t) = sum(sum(tss(:,1:Jw).*population(:,1:Jw))) / ...
        sum(sum(population(:,Jw+1:end)));
    TaxInc(t) = WMean(tax(:), population(:));
    Tax(t) = TaxInc(t) + TaxC(t);
    G(t) = Tax(t);
    % accidental bequest is the saving by dead agents
    Tr(t) = sum(sum(k(:,2:J).*(1-survival(:,1:J-1)).*population(:,1:J-1))) / ...
        sum(population(:));
    Y(t) = K(t)^Gamma * L(t)^(1-Gamma);
    GYRatio(t) = G(t) / Y(t);
    KYRatio(t) = K(t) / Y(t);
    KLRatio(t) = K(t) / L(t);
    R(t) = Gamma * KLRatio(t)^(Gamma-1) - Delta;
    W(t) = (1-Gamma) * KLRatio(t)^Gamma;
    MeanHours(t) = WMean(reshape(n(:,1:Jw)+s(:,1:Jw),1,[]),reshape(population(:,1:Jw),1,[]));
end

% write J+1 capital and population to empty
% SimulateResult.k(:,end) = [];
% SimulateResult.population(:,end) = [];

SimulateResult = v2struct(C,TaxC,K,L,Tss,TaxInc,Tax,G,Tr,Y,GYRatio,KYRatio,KLRatio,R,W,MeanHours,vFull,kFull,hFull);
ExitFlag = 0;
end

function s = WMean(x, weight)
s = sum(x.*weight)./sum(weight);
end