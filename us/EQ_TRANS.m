function EqResult = EQ_TRANS(VfiSs1Result, SimulateSs0Result, Params)
Params.TTrans = 150;
TTrans = Params.TTrans;
Rt = Params.r * ones(1,TTrans);
Wt = Params.w * ones(1,TTrans);
Trt = Params.Tr * ones(1,TTrans);
Tsst = Params.Tss * ones(1,TTrans);

Err = 1;
while (Err>1e-5)
    display('Vfi trans...');
    tic;
    VfiTransResult = Params.VfiFun(Params, VfiSs1Result, Rt, Wt, Trt, Tsst);
    toc;
    display('Simulate trans...');
    tic;
    SimulateTransResult = SIMULATE_TRANS(VfiTransResult, SimulateSs0Result, Params, Rt, Wt, Trt, Tsst);
    toc;
    NewRt = SimulateTransResult.R;
    NewWt = SimulateTransResult.W;
    NewTrt = SimulateTransResult.Tr;
    NewTsst = SimulateTransResult.Tss;
    
    Err = max( ...
        [abs(NewRt-Rt) abs(NewWt-Wt) abs(NewTrt-Trt) abs(NewTsst-Tsst)]);
    display(Err);
    
    Lambda = 0.2;
    Rt = (1-Lambda)*Rt + Lambda*NewRt;
    Wt = (1-Lambda)*Wt + Lambda*NewWt;
    Trt = (1-Lambda)*Trt + Lambda*NewTrt;
    Tsst = (1-Lambda)*Tsst + Lambda*NewTsst;
    
    save('SimulateTransResultBase', 'SimulateTransResult');
end
EqResult = SimulateTransResult;
end