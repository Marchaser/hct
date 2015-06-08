function EqResult = PARTIAL(Params)
Params = COMMON(Params);
tic;
[VfiResult, VfiExitFlag] = Params.VfiFun(Params);
toc;
tic;
[SimulateResult, SimulateExitFlag] = Params.SimulateFun(VfiResult,Params);
toc;
tic;
[AggregateResult, AggregateExitFlag] = AGGREGATE(SimulateResult,Params);
toc;

% PLOT(AggregateResult, Params);
EqResult = v2struct(Params, VfiResult, SimulateResult, AggregateResult);