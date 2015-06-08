function [EqResult, EqExitFlag] = BEQ(Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Solve balanced equilibrium. Tau0N is used to balance government budget.

% Params.AggregateDisplay = false;
LastX = [];
EqObjectsModel = [];

    function EqMetric = ComputeEqMetric(X)
        display(X);
        
        NewParams = Params;
        NewParams.Tau0N = X(1);
        NewParams.KLRatio = X(2);
        NewParams.Tss = X(3);
        NewParams.Tr = X(4);
        
        NewParams.r = Params.Gamma * NewParams.KLRatio^(Params.Gamma-1) - Params.Delta;
        NewParams.w = (1-Params.Gamma) * NewParams.KLRatio^Params.Gamma;
        
%         VfiTic = tic;
        [VfiResult, VfiExitFlag] = VFI_FACTOR_TAX(NewParams);
        if VfiExitFlag<0
            EqMetric = nan;
            return;
        end
%         toc(VfiTic);
%         SimulateTic = tic;
        [SimulateResult, SimulateExitFlag] = SIMULATE_FACTOR_TAX(VfiResult, NewParams);
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
        LastX = X;
        
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

options = optimoptions('fsolve','Display','iter','DiffMinChange',1e-3,'TolFun',1e-2,'OutputFcn',@outfun);
problem.objective = @ComputeEqMetric;
problem.x0 = Params.BEqX0;
problem.solver = 'fsolve';
problem.options = options;
fsolve(problem);

% return
EqResult = v2struct(VfiResult, SimulateResult, AggregateResult, EqObjectsModel, LastX);
EqExitFlag = 0;
end