Params = SETUP;
Params = MOMENTS(Params);
Params = COMMON(Params);
tic;
[VfiResult, VfiExitFlag] = VFI_FACTOR_TAX(Params);
toc;
tic;
[SimulateResult, SimulateExitFlag] = SIMULATE_FACTOR_TAX(VfiResult,Params);
toc;
tic;
[AggregateResult, AggregateExitFlag] = AGGREGATE(SimulateResult,Params);
toc;
save('AggregateResult', 'AggregateResult');

load('AggregateResult');
Params = MOMENTS(Params);
PLOT(AggregateResult, Params);