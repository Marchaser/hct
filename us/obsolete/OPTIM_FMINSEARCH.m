function OPTIM1(Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Solve optimal tax code using direct search. Tau2N is used to balance government budget.

EqResult = [];
WarmUpEqX0 = Params.BEqX0;

    function Welfare = ComputeWelfare(Tau)
        if any(Tau<Params.TauLb | Tau>Params.TauUb)
            Welfare = -inf;
            return;
        end
        
        NewParams = Params;
        % overwrite moments to params
        NewParams.Tau0K = Tau(1);
        NewParams.Tau1K = Tau(2);
        NewParams.Tau2K = Tau(3);
        NewParams.Tau1N = Tau(4);
        NewParams.Tau2N = Tau(5);
        
        % overwrite warm up initial guess
        NewParams.BEqX0 = WarmUpEqX0;
        
        fprintf('Current parameters:\n');
        display(Tau);
        display(WarmUpEqX0);
        
        EqTic = tic;
        [EqResult, EqExitFlag] = BEQ(NewParams);
        fprintf('Time for EQ:\n');
        toc(EqTic);
        
        Welfare = EqResult.AggregateResult.Welfare;
        display(Welfare);
        
        WarmUpEqX0 = EqResult.LastX;
    end

    function [f,g] = ComputeWelfareWrapper(x)
        f = -ComputeWelfare(x);
        if nargout>1
            [f,g] = AutoDiff(@(xx) -ComputeWelfare(xx), x, 1e-4, 'central');
        end
    end

x0 = [Params.TauX0];
lb = [Params.TauLb];
ub = [Params.TauUb];

% options = optimset('Algorithm', 'active-set', 'Display', 'iter', 'GradObj', 'on', 'TolFun', 1e-6);
% OptimTau = knitromatlab(@ComputeWelfareWrapper, x0, [], [], [], [], lb, ub, [], [], options);
OptimTau = fminsearch(@(x) -ComputeWelfare(x), x0);
end