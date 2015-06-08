function OPTIM1(Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Solve optimal tax code using fsolve with BEQ as nonlinear constraint.

Params.Tau1K = 0;

XHistory = [];
WelfareHistory = [];
EqMetricHistory = [];
    function Welfare = ComputeWelfare(X)
        % look up from history
        for i=1:size(XHistory,1)
            if isequal(X, XHistory(i,:))
                display('X found in history');
                Welfare = WelfareHistory(i);
                return;
            end
        end
        
        display('X not found in history');
        
        NewParams = Params;
        % overwrite moments to params
        NewParams.Tau0K = X(1);
        NewParams.Tau1N = X(2);
        NewParams.Tau2N = X(3);
        NewParams.Tau0N = X(4);
        NewParams.KLRatio = X(5);
        NewParams.Tss = X(6);
        NewParams.Tr = X(7);
        NewParams.r = Params.Gamma * NewParams.KLRatio^(Params.Gamma-1) - Params.Delta;
        NewParams.w = (1-Params.Gamma) * NewParams.KLRatio^Params.Gamma;
        
        fprintf('Current parameters:\n');
        display(X);
        
        [VfiResult, VfiExitFlag] = Params.VfiFun(NewParams);
        [SimulateResult, SimulateExitFlag] = Params.SimulateFun(VfiResult, NewParams);
        AggregateResult = AGGREGATE(SimulateResult, NewParams);
        
        Welfare = AggregateResult.Welfare;
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
        EqMetric = EqObjectsTarget - EqObjectsModel;
        
        % store in history
        XHistory = [XHistory;X];
        WelfareHistory = [WelfareHistory;Welfare];
        EqMetricHistory = [EqMetricHistory;EqMetric];
    end

    function EqMetric = ComputeEqMetric(X)        
        % look up from history
        for i=1:size(XHistory,1)
            if isequal(X, XHistory(i,:))
                display('X found in history');
                EqMetric = EqMetricHistory(i,:);
                return;
            end
        end
        
        display('X not found in history');
        
        NewParams = Params;
        % overwrite moments to params
        NewParams.Tau0K = X(1);
        NewParams.Tau1N = X(2);
        NewParams.Tau2N = X(3);
        NewParams.Tau0N = X(4);
        NewParams.KLRatio = X(5);
        NewParams.Tss = X(6);
        NewParams.Tr = X(7);
        NewParams.r = Params.Gamma * NewParams.KLRatio^(Params.Gamma-1) - Params.Delta;
        NewParams.w = (1-Params.Gamma) * NewParams.KLRatio^Params.Gamma;
        
        fprintf('Current parameters:\n');
        display(X);
        
        [VfiResult, VfiExitFlag] = Params.VfiFun(NewParams);
        [SimulateResult, SimulateExitFlag] = Params.SimulateFun(VfiResult, NewParams);
        AggregateResult = AGGREGATE(SimulateResult, NewParams);
        
        Welfare = AggregateResult.Welfare;
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
        EqMetric = EqObjectsTarget - EqObjectsModel;
        
        % store in history
        XHistory = [XHistory;X];
        WelfareHistory = [WelfareHistory;Welfare];
        EqMetricHistory = [EqMetricHistory;EqMetric];
    end

    function [f,g] = ComputeWelfareWrapper(x)
        f = -ComputeWelfare(x);
        if nargout>1
            [f,g] = AutoDiff(@(xx) -ComputeWelfare(xx), x, 1e-3, 'central');
        end
    end

    function [Empty1, Ceq, Empty2, CeqJac]  = ComputeEqMetricWrapper(x)
        Empty1 = [];
        Ceq = ComputeEqMetric(x);
        if nargout>2
            Empty2 = [];
            [Ceq,CeqJac] = AutoDiff(@ComputeEqMetric, x, 1e-3, 'central');
            CeqJac = CeqJac';
        end
    end

    function [stop,options,optchanged] = outfun(x,optimValues,state)
        stop = false;
        options = [];
        optchanged = [];
        fprintf('Optim iteration\n');
        display(x');
    end

x0 = [Params.TauX0([1 4 5]) Params.BEqX0];
lb = [Params.TauLb([1 4 5]) Params.BEqLb];
ub = [Params.TauUb([1 4 5]) Params.BEqUb];

JPattern = sparse(ones(4, 7));
ExtendedFeatures.JacobPattern = JPattern;
options = optimset('Algorithm', 'interior-point', 'Display', 'iter', 'GradObj', 'on', 'GradConstr', 'on', 'TolFun', 1e-6, 'TolCon', 1e-6);

OptimTau = knitromatlab(@ComputeWelfareWrapper, x0, [], [], [], [], lb, ub, @ComputeEqMetricWrapper, [], options);
end
