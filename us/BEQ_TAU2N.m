function [EqResult, EqExitFlag] = BEQ(Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Solve balanced equilibrium. Tau2N is used to balance government budget.
% Use lsqnonlin to enforce boundary

Params.BEqX0(1) = [];
Params.BEqLb(1) = [];
Params.BEqUb(1) = [];
Tau2NAtEq = Params.Tau2N;

% Params.AggregateDisplay = false;
LastX = [];
EqObjectsModel = [];
EqMetricHistory = [];
XHistory = [];

    function EqMetric = ComputeEqMetric(X)
        % serach in history
        for i=1:size(XHistory,1)
            if isequal(X', XHistory(i,:));
                EqMetric = EqMetricHistory(i,:);
                return;
            end
        end
        display(X');
        
        NewParams = Params;
        NewParams.KLRatio = X(1);
        NewParams.Tss = X(2);
        NewParams.Tr = X(3);
        
        NewParams.r = Params.Gamma * NewParams.KLRatio^(Params.Gamma-1) - Params.Delta;
        NewParams.w = (1-Params.Gamma) * NewParams.KLRatio^Params.Gamma;
        
        AggregateResult = [];
        
        function GDiff = ComputeEqMetricVaryTau2N(x)
            % solve a nested problem using Tau2N to target G
            NewNewParams = NewParams;
            NewNewParams.Tau2N = x;
            display(x);
            
            
            %         VfiTic = tic;
            [VfiResult, VfiExitFlag] = Params.VfiFun(NewNewParams);
            if VfiExitFlag<0
                GDiff = nan;
                return;
            end
            %         toc(VfiTic);
            %         SimulateTic = tic;
            [SimulateResult, SimulateExitFlag] = Params.SimulateFun(VfiResult, NewNewParams);
            if SimulateExitFlag<0
                GDiff = nan;
                return;
            end
            %         toc(SimulateTic);
            AggregateResult = AGGREGATE(SimulateResult, NewNewParams);
            GDiff = NewNewParams.GBench - AggregateResult.G;
        end
        
        Tau2NAtEq = fzero(@ComputeEqMetricVaryTau2N, [0 1e9], optimset('TolX',1e-3));
        
        EqObjectsTarget = [
            NewParams.KLRatio
            NewParams.Tss
            NewParams.Tr
            ]';
        
        EqObjectsModel = [
            AggregateResult.KLRatio
            AggregateResult.Tss
            AggregateResult.Tr
            ]';
        
        display(EqObjectsTarget);
        display(EqObjectsModel);
        
        EqMetric = EqObjectsTarget - EqObjectsModel;
        
        XHistory = [X' ; XHistory];
        EqMetricHistory = [EqMetric ; EqMetricHistory];
    end

    function [stop,options,optchanged] = outfun(x,optimValues,state)
        stop = false;
        options = [];
        optchanged = [];
        fprintf('Eq iteration\n');
    end

%{
options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','Display','iter','DiffMinChange',1e-3,'TolFun',1e-3,'OutputFcn',@outfun);
problem.objective = @ComputeEqMetric;
problem.x0 = Params.BEqX0;
problem.lb = Params.BEqLb;
problem.ub = Params.BEqUb;
problem.solver = 'lsqnonlin';
problem.options = options;
LastX = lsqnonlin(problem);
display('BEQ solved');
%}

%{
% Check feasibility
MaxTaxBeqX = Params.BEqX0';
MaxTaxBeqX(1) = 1e5;
MaxEqMetric = ComputeEqMetric(MaxTaxBeqX);
MaxTax = MaxEqMetric(1);
if MaxTax>0
    % not feasible, return NaN in welfare
    EqResult = [];
    EqExitFlag = -1;
    return;
end
display('Laffer Curve Feasible');

% First solve Tau2N to a low precision
    function GMetric = ComputeEqMetricVaryTau2N(Tau2N)
        EqMetric = ComputeEqMetric([Tau2N Params.BEqX0(2:end)]');
        GMetric = EqMetric(1);
    end
options = [];
Tau2NApprox = CoDoSol(Params.BEqX0(1), @ComputeEqMetricVaryTau2N, 0, 1e5, [1e-3, 0], options)';
display('Tau2NApprox solved');

% Find balanced equilibrium
Params.BEqX0(1) = Tau2NApprox;
%}

options = [];
LastX = CoDoSol(Params.BEqX0', @(x) ComputeEqMetric(x)', Params.BEqLb', Params.BEqUb', [1e-5, 0], options)';
display(LastX);
display('BEQ solved');

% return
EqResult = v2struct(VfiResult, SimulateResult, AggregateResult, EqObjectsModel, LastX);
EqExitFlag = 0;
end