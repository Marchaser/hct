function CALIBRATE_NO_TYPE(Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Calibrate model. The idea is that static moments (KYRatio, Mean
% Hours etc.) are matched along with finding equilibrium implemented in EQ.
% Then human capital process is calibrated to match the age profiles.
% Use lsqnonlin

Params.LogH1Sigma = 1e-12;
Params.LogASigma = 1e-12;
Params.CaliX0 = [Params.Alpha Params.LogAMu];
Params.CaliLb = [0 -2];
Params.CaliUb = [1 -1];

EqResult = [];
WarmUpEqX = Params.EqX0;
LastEqX = WarmUpEqX;
MinMetricFound = 999;

    function MomentsMetric = ComputeMomentsMetric(X)
        NewParams = Params;
        NewParams.Alpha = X(1);
        NewParams.LogAMu = X(2);
        
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
        
        MomentsMetric = EqResult.AggregateResult.MomentsModel(1:40) - Params.MomentsData(1:40);
        
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

options = optimoptions('lsqnonlin','DiffMinChange',1e-3,'Display','iter');
lsqnonlin(@ComputeMomentsMetric, Params.CaliX0, Params.CaliLb, Params.CaliUb, options);
end

