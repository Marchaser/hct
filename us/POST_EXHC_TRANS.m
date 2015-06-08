load('SimulateTransResultExhc');
Params = SETUP;
Params = OVERWRITE(Params, 'Benchmark', 'Exhc');
Params = COMMON(Params);
Params.ControlSpec = 'Exhc';
Params.TaxSpec = 'Factor';
VfiSs1Result = VFI_EXHC_FACTOR(Params);
Params.TTrans = 100;
Rt = SimulateTransResult.R;
Wt = SimulateTransResult.W;
Trt = SimulateTransResult.Tr;
Tsst = SimulateTransResult.Tss;
% VfiTransResult = VFI_EXHC_FACTOR_TRANS(Params, VfiSs1Result, Rt, Wt, Trt, Tsst);
% save('VfiTransResultExhc','VfiTransResult');
load('VfiTransResultExhc');
Params.TTrans = 1;
SimulateTransFullResult = SIMULATE_TRANS_FULL(VfiTransResult,SimulateSs0Result,Params,Rt,Wt,Trt,Tsst);
save('SimulateTransFullResultExhc','SimulateTransFullResult');

Params = SETUP;
Params = OVERWRITE(Params, 'Benchmark');
Params = COMMON(Params);
VfiSs0Result = VFI(Params);
SimulateSs0Result = SIMULATE(VfiSs0Result, Params);
AggregateSs0Result = AGGREGATE(SimulateSs0Result, Params);
SimulateSs0FullResult = SIMULATE_FULL(VfiSs0Result, Params);
save('SimulateSs0FullResult', 'SimulateSs0FullResult');

vDiff = SimulateTransFullResult.vFull(:,:,1) - SimulateSs0FullResult.v;
[vDiff(:,2) SimulateSs0FullResult.k(:,2)];
SupportByJ = mean(vDiff>0 .* SimulateSs0FullResult.population);
TotalSupport = sum(double(vDiff(:)>0) .* SimulateSs0FullResult.population(:)) / sum(SimulateSs0FullResult.population(:));
OldDie = sum(sum(SimulateSs0FullResult.population(:,46:end))) / sum(sum(SimulateSs0FullResult.population(:,1:end)));
% support by a and j
for j=1:Params.Jw
    SupportByJA1(j) = mean(vDiff(Params.AIdx==1,j)>0 .* SimulateSs0FullResult.population(Params.AIdx==1,j));
    SupportByJA2(j) = mean(vDiff(Params.AIdx==2,j)>0 .* SimulateSs0FullResult.population(Params.AIdx==2,j));
end
figure;
hold on;
plot(SupportByJA1,'k');
plot(SupportByJA2,'b--');
hold off;

figure;
hold on;
plot(SimulateTransResult.Tss);
hold off