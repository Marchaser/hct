clear;

Params = SETUP;
Params = OVERWRITE(Params, 'Benchmark', 'Base');
Params.ControlSpec = 'Base';
Params.TaxSpec = 'Factor';
Params = COMMON(Params);
% VfiNoRiskResult = VFI_NO_RISK(Params);
% SimulateNoRiskResult = SIMULATE(VfiNoRiskResult, Params);
% AggregateNoRiskResult = AGGREGATE(SimulateNoRiskResult, Params);

% Params.Tau0K = Params.Tau0K + 0.01;
% Params.Tau1N = Params.Tau1N;
VfiResult = VFI_FACTOR(Params);
SimulateResult = SIMULATE(VfiResult, Params);
AggregateResult = AGGREGATE(SimulateResult, Params);