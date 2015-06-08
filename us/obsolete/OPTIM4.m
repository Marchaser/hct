function OPTIM4(Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Solve optimal tax code using fmincon. Tau2N is used to balance government budget.

EqResult = [];
WarmUpEqX0 = Params.BEqX0;

    function Welfare = ComputeWelfare(Tau)
        NewParams = Params;
        % overwrite moments to params
        NewParams.Tau0K = Tau(1);
        NewParams.Tau1N = Tau(2);
        NewParams.Tau2N = Tau(3);
        
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
        
        % display something
        display(EqResult.LastX);
        display(EqResult.EqObjectsModel);
        
        WarmUpEqX0 = EqResult.LastX;
    end

    function [stop,options,optchanged] = outfun(x,optimValues,state)
        stop = false;
        options = [];
        optchanged = [];
        fprintf('Optim iteration\n');
    end

options = optimoptions('fmincon','Display','Iter','DiffMinChange',1e-2,'Algorithm','sqp','OutputFcn',@outfun);
fmincon(@(x) -ComputeWelfare(x), Params.TauReducedX0, [], [], [], [], Params.TauReducedLb, Params.TauReducedUb, [], options);
end