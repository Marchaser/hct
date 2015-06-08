Params = SETUP;
Params = OVERWRITE(Params, 'Benchmark');
Params = COMMON(Params);
VfiSs0Result = VFI(Params);
SimulateSs0Result = SIMULATE(VfiSs0Result, Params);
AggregateSs0Result = AGGREGATE(SimulateSs0Result, Params);

Params = SETUP;
Params = OVERWRITE(Params, 'Benchmark', 'Base');
Params = COMMON(Params);
VfiSs1Result = VFI_FACTOR(Params);
Params.VfiFun = @VFI_FACTOR_TRANS;
EQ_TRANS(VfiSs1Result, SimulateSs0Result, Params);