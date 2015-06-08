function [EqResult, ExitFlag] = EQ(Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Finding equilibrium using fsolve. Several static moments (KYRatio, Mean
% Hours etc.) are matched.

% Params.AggregateDisplay = false;
LastX = [];
EqObjectsModel = [];

    function EqMetric = ComputeEqMetric(X)
        LastX = X';
        
        NewParams = Params;
        NewParams.Beta = X(1);
        NewParams.Chi = X(2);
        NewParams.Tau2Gs = X(3);
        NewParams.Tss = X(4);
        NewParams.Tr = X(5);
%         NewParams.KLRatio = X(6);
%         NewParams.r = NewParams.Gamma * NewParams.KLRatio^(NewParams.Gamma-1) - NewParams.Delta;
%         NewParams.w = (1-NewParams.Gamma) * NewParams.KLRatio^NewParams.Gamma;
        
%         VfiTic = tic;
        VfiResult = VFI(NewParams);
%         toc(VfiTic);
%         SimulateTic = tic;
        SimulateResult = SIMULATE(VfiResult, NewParams);
%         toc(SimulateTic);
        AggregateResult = AGGREGATE(SimulateResult, NewParams);
        
        EqObjectsTarget = [
            NewParams.KYRatioTarget
            NewParams.MeanHoursTarget
            NewParams.GYRatioTarget
            NewParams.Tss
            NewParams.Tr
            ]';
        EqObjectsModel = [
            AggregateResult.KYRatio
            AggregateResult.MeanHours
            AggregateResult.GYRatio
            AggregateResult.Tss
            AggregateResult.Tr
            ]';
        
        EqMetric = EqObjectsTarget - EqObjectsModel;
    end

%{
options = optimoptions('fsolve','Display','iter','DiffMinChange',1e-3,'TolFun',1e-2,'TolX',1e-6);
problem.objective = @ComputeEqMetric;
problem.x0 = Params.EqX0;
problem.solver = 'fsolve';
problem.options = options;
fsolve(problem);
%}

options = [];
LastX = CoDoSol(Params.EqX0', @(x) ComputeEqMetric(x)', Params.EqLb', Params.EqUb', [1e-4, 0], options)';
display(LastX);

% return
EqResult = v2struct(VfiResult, SimulateResult, AggregateResult, EqObjectsModel, LastX);
ExitFlag = 0;
end
