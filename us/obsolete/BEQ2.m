function [EqResult, EqExitFlag] = BEQ(Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Solve balanced equilibrium. Tau0N is used to balance government budget.
% Use lsqnonlin to enforce boundary

% Params.AggregateDisplay = false;
LastX = [];
EqObjectsModel = [];

    function EqMetric = ComputeEqMetric(X)
        display(X');
        
        NewParams = Params;
        NewParams.Tau0N = X(1);
        NewParams.KLRatio = X(2);
        NewParams.Tss = X(3);
        NewParams.Tr = X(4);
        
        NewParams.r = Params.Gamma * NewParams.KLRatio^(Params.Gamma-1) - Params.Delta;
        NewParams.w = (1-Params.Gamma) * NewParams.KLRatio^Params.Gamma;
        
%         VfiTic = tic;
        [VfiResult, VfiExitFlag] = Params.VfiFun(NewParams);
        if VfiExitFlag<0
            EqMetric = nan;
            return;
        end
%         toc(VfiTic);
%         SimulateTic = tic;
        [SimulateResult, SimulateExitFlag] = Params.SimulateFun(VfiResult, NewParams);
        if SimulateExitFlag<0
            EqMetric = nan;
            return;
        end
%         toc(SimulateTic);
        AggregateResult = AGGREGATE(SimulateResult, NewParams);
        
        EqObjectsTarget = [
            NewParams.GBench
            NewParams.KLRatio
            NewParams.Tss
            NewParams.Tr
            ]';
        
        EqObjectsModel = [
            AggregateResult.G
            AggregateResult.KLRatio
            AggregateResult.Tss
            AggregateResult.Tr
            ]';
        
        display(EqObjectsTarget);
        display(EqObjectsModel);
        
        EqMetric = EqObjectsTarget - EqObjectsModel;
    end

    function [stop,options,optchanged] = outfun(x,optimValues,state)
        stop = false;
        options = [];
        optchanged = [];
        fprintf('Eq iteration\n');
    end

%{
options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','iter','DiffMinChange',1e-3,'TolFun',1e-3,'OutputFcn',@outfun);
problem.objective = @ComputeEqMetric;
problem.x0 = Params.BEqX0;
problem.lb = Params.BEqLb;
problem.ub = Params.BEqUb;
problem.solver = 'lsqnonlin';
problem.options = options;
lsqnonlin(problem);
%}

options = [];
LastX = CoDoSol(Params.BEqX0', @(x) ComputeEqMetric(x)', Params.BEqLb', Params.BEqUb', [1e-5, 0], options)';
display(LastX);

% return
EqResult = v2struct(VfiResult, SimulateResult, AggregateResult, EqObjectsModel, LastX);
EqExitFlag = 0;
end