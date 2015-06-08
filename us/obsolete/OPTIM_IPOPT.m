function OPTIM1(Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Solve optimal tax code using fsolve with BEQ as nonlinear constraint.
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
        NewParams.Tau1K = X(2);
        NewParams.Tau2K = X(3);
        NewParams.Tau1N = X(4);
        NewParams.Tau2N = X(5);
        NewParams.Tau0N = X(6);
        NewParams.KLRatio = X(7);
        NewParams.Tss = X(8);
        NewParams.Tr = X(9);
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
        NewParams.Tau1K = X(2);
        NewParams.Tau2K = X(3);
        NewParams.Tau1N = X(4);
        NewParams.Tau2N = X(5);
        NewParams.Tau0N = X(6);
        NewParams.KLRatio = X(7);
        NewParams.Tss = X(8);
        NewParams.Tr = X(9);
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

    function g = ComputeWelfareGradient(x)
        [~,g] = AutoDiff(@(xx) -ComputeWelfare(xx), x, 1e-4, 'central');
    end

    function j = ComputeEqMetricJacobian(x)
        [~,j] = AutoDiff(@ComputeEqMetric, x, 1e-4, 'central');
        j = sparse(j);
        display('End of Jacobian');
    end

funcs.objective = @(x) -ComputeWelfare(x);
funcs.gradient = @ComputeWelfareGradient;
funcs.constraints = @ComputeEqMetric;
funcs.jacobian = @ComputeEqMetricJacobian;
funcs.jacobianstructure = @(x) sparse(ones(4,9));

x0 = [Params.TauX0 Params.BEqX0];
options.lb = [Params.TauLb Params.BEqLb];
options.ub = [Params.TauUb Params.BEqUb];
options.cl = zeros(1,4);
options.cu = zeros(1,4);
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.limited_memory_max_history = 1e4;

[OptimTau, info] = ipopt(x0, funcs, options);
end