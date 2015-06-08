clear;

%{
Params = SETUP;
Params = OVERWRITE(Params, 'Benchmark');
Params = COMMON(Params);
VfiSs0Result = VFI(Params);
SimulateSs0Result = SIMULATE(VfiSs0Result, Params);
AggregateSs0Result = AGGREGATE(SimulateSs0Result, Params);

%{
Params = SETUP;
Params = OVERWRITE(Params, 'Benchmark', 'Exhc');
Params.ControlSpec = 'Exhc';
Params.TaxSpec = 'Factor';
Params.VfiFun = @VFI_EXHC_FACTOR;
ExhcAtOptim = PARTIAL(Params);
%}

Params = SETUP;
Params = OVERWRITE(Params, 'Benchmark', 'Base');
Params = COMMON(Params);
VfiSs1Result = VFI_FACTOR(Params);

load('VfiTransResult_2');
load('SimulateTransResult');
Params.TTrans = 100;
Rt = SimulateTransResult.R;
Wt = SimulateTransResult.W;
Trt = SimulateTransResult.Tr;
Tsst = SimulateTransResult.Tss;
VfiTransResult = VFI_FACTOR_TRANS(Params, VfiSs1Result, Rt, Wt, Trt, Tsst);
% save('VfiTransResult','VfiTransResult');

Params.TTrans = 1;
SimulateSs0FullResult = SIMULATE_FULL(VfiSs0Result, Params);
SimulateTransFullResult = SIMULATE_TRANS_FULL(VfiTransResult,SimulateSs0Result,Params,Rt,Wt,Trt,Tsst);

clear;
Params = SETUP;
Params = OVERWRITE(Params, 'Benchmark', 'Exhc');
Params = COMMON(Params);
Params.ControlSpec = 'Exhc';
Params.TaxSpec = 'Factor';
VfiSs1Result = VFI_EXHC_FACTOR(Params);
SimulateSs1Result = SIMULATE(VfiSs1Result,Params);
AggregateSs1Result = AGGREGATE(SimulateSs1Result,Params);
Params.TTrans = 1;
VfiTransResult = VFI_EXHC_FACTOR_TRANS(Params,VfiSs1Result,Params.r,Params.w,Params.Tr,Params.Tss);
SimulateTransResult = SIMULATE_TRANS(VfiTransResult,SimulateSs0Result,Params,Params.r,Params.w,Params.Tr,Params.Tss);
SimulateTransResult = SIMULATE_TRANS(VfiTransResult,SimulateSs1Result,Params,Params.r,Params.w,Params.Tr,Params.Tss);
% EQ_TRANS(VfiSs1Result, SimulateSs0Result, Params);
%}

Params = SETUP;
Params = OVERWRITE(Params, 'Benchmark');
Params = COMMON(Params);
VfiSs0Result = VFI(Params);
SimulateSs0Result = SIMULATE(VfiSs0Result, Params);
AggregateSs0Result = AGGREGATE(SimulateSs0Result, Params);

Params = SETUP;
Params = OVERWRITE(Params, 'Benchmark', 'Exhc');
Params = COMMON(Params);
Params.ControlSpec = 'Exhc';
Params.TaxSpec = 'Factor';
VfiSs1Result = VFI_EXHC_FACTOR(Params);
SimulateSs1Result = SIMULATE(VfiSs1Result, Params);
AggregateSs1Result = AGGREGATE(SimulateSs1Result,Params);

Params.TTrans = 1;
VfiTransResult = VFI_EXHC_FACTOR_TRANS(Params,VfiSs1Result,Params.r,Params.w,Params.Tr,Params.Tss);
SimulateTransResult = SIMULATE_TRANS(VfiTransResult,SimulateSs1Result,Params,Params.r,Params.w,Params.Tr,Params.Tss);

Params.VfiFun = @VFI_EXHC_FACTOR_TRANS;
EqResult = EQ_TRANS(VfiSs1Result, SimulateSs0Result, Params);


Params.TTrans = 100;
Rt = SimulateTransResult.R;
Wt = SimulateTransResult.W;
Trt = SimulateTransResult.Tr;
Tsst = SimulateTransResult.Tss;
VfiTransResult = VFI_EXHC_FACTOR_TRANS(Params, VfiSs1Result, Rt, Wt, Trt, Tsst);
% save('VfiTransResult_exhc','VfiTransResult');

Params.TTrans = 1;
SimulateTransFullResult = SIMULATE_TRANS_FULL(VfiTransResult,SimulateSs0Result,Params,Rt,Wt,Trt,Tsst);

Params = SETUP;
Params = OVERWRITE(Params, 'Benchmark');
Params = COMMON(Params);
SimulateSs0FullResult = SIMULATE_FULL(VfiSs0Result, Params);

% 
% Rt = Params.r * ones(1,5);
% Wt = Params.w * ones(1,5);
% Trt = Params.Tr * ones(1,5);
% Tsst = Params.Tss * ones(1,5);
% 
% 
% VfiTransResult = VFI_FACTOR_TRANS(Params, VfiSs1Result, Rt, Wt, Trt, Tsst);
% SimulateTransResult = SIMULATE_TRANS(VfiTransResult, SimulateSs0Result, Params, Rt, Wt, Trt, Tsst);