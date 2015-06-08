function OPTIM2(Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Solve optimal tax code using direct search. Tau2N is used to balance government budget.

EqResult = [];
WarmUpEqX0 = Params.BEqX0;
LastEqX0 = WarmUpEqX0;

MaxWelfareFound = -999;

    function Welfare = ComputeWelfare(Tau)
        % scale back
        Tau = Tau .* Params.TauScale;
        
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
        
        % display something
        display(EqResult.LastX);
        display(EqResult.EqObjectsModel);
        
        if (Welfare > MaxWelfareFound)
            LastEqX0 = EqResult.LastX;
            MaxWelfareFound = Welfare;
        end
    end

    function [stop,options,optchanged] = outfun(x,optimValues,state)
        stop = false;
        options = [];
        optchanged = [];
        WarmUpEqX0 = LastEqX0;
        fprintf('Optim iteration\n');
%         switch state
%             case 'iter'
%                 WarmUpEqObjects = LastEqObjects;
%             otherwise
%         end
    end

rng(0823);
options = psoptimset('MaxMeshSize', 1, 'InitialMeshSize', 1, 'Display', 'iter', 'PollingOrder', 'Consecutive', 'ScaleMesh', 'Off','CompletePoll','Off','OutputFcns',@outfun);
patternsearch(@(x) -ComputeWelfare(x), Params.TauX0./Params.TauScale, [], [], [], [], Params.TauLb./Params.TauScale, Params.TauUb./Params.TauScale, [], options);

end