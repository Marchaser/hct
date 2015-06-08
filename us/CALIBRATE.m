function CALIBRATE(Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Calibrate model. The idea is that static moments (KYRatio, Mean
% Hours etc.) are matched along with finding equilibrium implemented in EQ.
% Then human capital process is calibrated to match the age profiles.
% Use lsqnonlin

EqResult = [];
WarmUpEqX = Params.EqX0;
LastEqX = WarmUpEqX;
MinMetricFound = 999;

switch (Params.Profile)
    case 'NoIdio'
        Params.CaliX0 = [Params.LogEpsilonMu];
        Params.CaliLb = [-1];
        Params.CaliUb = [1];
    case 'NoType'
        Params.CaliX0 = [Params.LogAMu Params.LogH1Mu];
        Params.CaliLb = [-2 -1];
        Params.CaliUb = [-1 1];
end

    function MomentsMetric = ComputeMomentsMetric(X)
        NewParams = Params;
        switch(Params.Profile)
            case 'BenchMark'
                NewParams.Alpha = X(1);
                NewParams.LogH1Sigma = X(2);
                NewParams.LogAMu = X(3);
                NewParams.LogASigma = X(4);
                NewParams.HARho = X(5);
            case 'NoIdio'
                NewParams.LogEpsilonMu = X(1);
            case 'NoType'
                NewParams.LogAMu = X(1);
                NewParams.LogH1Mu = X(2);
        end
        
        % overwrite warm up equilibrium initial values
        NewParams.EqX0 = WarmUpEqX;
        
        % NOTE(wenlan): need to call COMMON here to regenerate random
        % numbers since distributions have changed!
        NewParams = COMMON(NewParams);
        
        fprintf('Current parameters:\n');
        display(X);
        display(WarmUpEqX);
        
        EqTic = tic;
        [EqResult, EqExitFlag] = EQ(NewParams);
        fprintf('Time for Eq:\n');
        toc(EqTic);
        
        switch(Params.Profile)
            case 'BenchMark'
                MomentsMetric = (EqResult.AggregateResult.MomentsModel - Params.MomentsData) ./ Params.MomentsData;
            case {'NoIdio', 'NoType'}
                MomentsMetric = EqResult.AggregateResult.MomentsModel(1:40) - Params.MomentsData(1:40);
        end
        
        % display something
%         display(MomentsMetric);
        display(EqResult.LastX);
        display(EqResult.EqObjectsModel);
        display([EqResult.AggregateResult.MomentsModel; Params.MomentsData]);
        
        WarmUpEqX = EqResult.LastX;
    end

    function [stop,options,optchanged] = outfun(x,optimValues,state)
        stop = false;
        options = [];
        optchanged = [];
        WarmUpEqX = LastEqX;
    end

options = optimoptions('lsqnonlin','FinDiffRelStep',1e-3,'FinDiffType','central','Display','iter');
lsqnonlin(@ComputeMomentsMetric, Params.CaliX0, Params.CaliLb, Params.CaliUb, options);
end

