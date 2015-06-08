function OPTIM3(Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Solve optimal tax code using fmincon. Tau2N is used to balance government budget.

EqResult = [];
WarmUpEqX0 = Params.BEqX0;
TauHistory = [];
WelfareHistory = [];

    function Welfare = ComputeWelfare(Tau)
%         Tau = Tau .* Params.TauScale;
        
        % find if Tau exists in history
        for i=1:size(TauHistory, 1)
            if isequal(Tau, TauHistory(i,:))
                Welfare = WelfareHistory(i);
                display('Tau found in history');
                return;
            end
        end
        
        display('Tau not found in history');

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
        if EqExitFlag<0
            Welfare = nan;
            return;
        end
        
        Welfare = EqResult.AggregateResult.Welfare;
        display(Welfare);
        
        % display something
        display(EqResult.LastX);
        display(EqResult.EqObjectsModel);
        
        WarmUpEqX0 = EqResult.LastX;
        
        % write in history
        TauHistory = [TauHistory; Tau];
        WelfareHistory = [WelfareHistory; Welfare];
    end

    function [f,g] = ComputeWelfareWrapper(x)
        if nargout==1
            f = -ComputeWelfare(x);
        elseif nargout==2
            [f,g] = AutoDiff(@(xx) -ComputeWelfare(xx), x, 1e-5);
        end
    end

    function [stop,options,optchanged] = outfun(x,optimValues,state)
        stop = false;
        options = [];
        optchanged = [];
        fprintf('Optim iteration\n');
        display(x);
    end

options = optimoptions('fmincon','Display','Iter','Algorithm','sqp','OutputFcn',@outfun, 'GradObj', 'on');
fmincon(@ComputeWelfareWrapper, Params.TauX0, [], [], [], [], Params.TauLb, Params.TauUb, [], options);
end