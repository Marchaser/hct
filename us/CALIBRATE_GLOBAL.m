function CALIBRATE_GLOBAL(Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Calibrate model. The idea is that static moments (KYRatio, Mean
% Hours etc.) are matched along with finding equilibrium implemented in EQ.
% Then human capital process is calibrated to match the age profiles.

EqResult = [];
WarmUpEqX = Params.EqX0;
LastEqX = WarmUpEqX;
MinMetricFound = 999;

    function MomentsMetric = ComputeMomentsMetric(X)
        % scale back
        X = X .* Params.CaliScale;
        % display parameters
        
        Newparams = Params;
        Newparams.Alpha = X(1);
        Newparams.LogH1Sigma = X(2);
        Newparams.LogAMu = X(3);
        Newparams.LogASigma = X(4);
        Newparams.H1ARho = X(5);
        
        % overwrite warm up equilibrium initial values
        Newparams.EqX0 = WarmUpEqX;
        
        % NOTE(wenlan): need to call COMMON here to regenerate random
        % numbers since distributions have changed!
        Newparams = COMMON(Newparams);
        
        fprintf('Current parameters:\n');
        display(X);
        display(WarmUpEqX);
        
        EqTic = tic;
        [EqResult, EqExitFlag] = EQ(Newparams);
        fprintf('Time for Eq:\n');
        toc(EqTic);
        
        MomentsMetric = sum((EqResult.AggregateResult.MomentsModel - Params.MomentsData).^2 .* Params.MomentsWeight);
        
        % display something
        display(MomentsMetric);
        display(EqResult.LastX);
        display(EqResult.EqObjectsModel);
        display([EqResult.AggregateResult.MomentsModel; Params.MomentsData]);

        if (MomentsMetric < MinMetricFound)
            LastEqX = EqResult.LastX;
            MinMetricFound = MomentsMetric;
        end
    end

    function [stop,options,optchanged] = outfun(x,optimValues,state)
        stop = false;
        options = [];
        optchanged = [];
        WarmUpEqX = LastEqX;
    end

options = psoptimset('MaxMeshSize', 1, 'InitialMeshSize', 1, 'Display', 'iter', 'PollingOrder', 'Consecutive', 'ScaleMesh', 'Off','CompletePoll','On','OutputFcns',@outfun);
patternsearch(@ComputeMomentsMetric, Params.CaliX0./Params.CaliScale, [], [], [], [], Params.CaliLb./Params.CaliScale, Params.CaliUb./Params.CaliScale, [], options);
end

