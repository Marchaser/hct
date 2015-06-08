% Task = 'EXHC';
% SpecIdio = true;
% SpecType = true;
% SpecBC = true;

format shortG;

if ~isequal(exist('Task','var'),1)
    Task = 'US';
end

Params = SETUP;
Params = MOMENTS(Params);

switch (Task)
    case 'EQ'
        Params = OVERWRITE(Params, 'Benchmark');
        Params = COMMON(Params);
        EQ(Params);
    case 'EQ_EXN'
        Params = OVERWRITE(Params, 'Benchmark', 'Exn');
        Params = COMMON(Params);
        EQ(Params);
    case 'US'
        Params = OVERWRITE(Params, 'Benchmark');
        BenchResult = PARTIAL(Params);
        save('BenchResult','BenchResult');
    case 'EXHC'
        Params = OVERWRITE(Params, 'Benchmark');
        Params.ControlSpec = 'Exhc';
        Params.VfiFun = @VFI_EXHC;
        Params = COMMON(Params);
        BenchResult = PARTIAL(Params);
    case 'EXN'
        Params = OVERWRITE(Params, 'Benchmark', 'Exn');
        Params = COMMON(Params);
        BenchResult = PARTIAL(Params);
    case 'CALIBRATE'
        Params = OVERWRITE(Params, 'Benchmark');
        CALIBRATE(Params);
    case 'CALIBRATE_EXN'
        Params = OVERWRITE(Params, 'Benchmark', 'Exn');
        CALIBRATE(Params);
    case 'OPTIM_BASE'
        Params = OVERWRITE(Params, 'Benchmark', 'Base');
        Params.ControlSpec = 'Base';
        Params.TaxSpec = 'Factor';
        Params = COMMON(Params);
        Params.VfiFun = @VFI_FACTOR;
        OPTIM1(Params);
    case 'OPTIM_EXN'
        Params = OVERWRITE(Params, 'Benchmark', 'Exn');
%         Params.ControlSpec = 'Base';
        Params.TaxSpec = 'Factor';
        Params = COMMON(Params);
        Params.VfiFun = @VFI_FACTOR;
        OPTIM1(Params);
    case 'OPTIM_EXHC'
        Params = OVERWRITE(Params, 'Benchmark', 'Exhc');
        Params.ControlSpec = 'Exhc';
        Params.TaxSpec = 'Factor';
        Params = COMMON(Params);
        Params.VfiFun = @VFI_EXHC_FACTOR;
        OPTIM1(Params);
    case 'OPTIM_EXALL'
        Params = OVERWRITE(Params, 'Benchmark', 'ExAll');
        Params.ControlSpec = 'Exhc';
        Params.TaxSpec = 'Factor';
        Params = COMMON(Params);
        Params.VfiFun = @VFI_EXHC_FACTOR;
        OPTIM1(Params);
    case 'BASE_AT_OPTIM'
        Params = OVERWRITE(Params, 'Benchmark', 'Base');
        Params.ControlSpec = 'Base';
        Params.TaxSpec = 'Factor';
        Params.VfiFun = @VFI_FACTOR;
        BaseAtOptim = PARTIAL(Params);
        save('BaseAtOptim','BaseAtOptim');
    case 'BASE_AT_SUBOPTIM'
        Params = OVERWRITE(Params, 'Benchmark', 'Exhc');
%         Params.BEqX0 = [2.4006e+08 4.3946        0.346     0.028036];
%         Params.Tau0N = 0.22246;
%         Params.Tau2N = 2.4006e+08;
%         Params.KLRatio = 4.3946;
%         Params.Tss = 0.346;
%         Params.Tr = 0.028036;
        Params.ControlSpec = 'Base';
        Params.TaxSpec = 'Factor';
        Params.VfiFun = @VFI_FACTOR;
        Params = COMMON(Params);
%         Params.Tau2N = 1e8;
        BenchResult = PARTIAL(Params);
        save('BaseAtSubOptim', 'BenchResult');
%         BEQ_TAU2N(Params);
    case 'EXHC_AT_OPTIM'
        Params = OVERWRITE(Params, 'Benchmark', 'Exhc');
        Params.ControlSpec = 'Exhc';
        Params.TaxSpec = 'Factor';
        Params.VfiFun = @VFI_EXHC_FACTOR;
        ExhcAtOptim = PARTIAL(Params);
        save('ExhcAtOptim','ExhcAtOptim');
    case 'EXN_AT_OPTIM'
        Params = OVERWRITE(Params, 'Benchmark', 'Exn');
        Params.ControlSpec = 'Base';
        Params.TaxSpec = 'Factor';
        Params.VfiFun = @VFI_FACTOR;
        BenchResult = PARTIAL(Params);
        save('ExnAtOptim','BenchResult');        
    case 'CALIBRATE_NO_BC'
        Params = OVERWRITE(Params, 'NoBc');
        CALIBRATE(Params);
    case 'NO_BC'
        Params = OVERWRITE(Params, 'NoBc');
        Params = COMMON(Params);
        BenchResult = PARTIAL(Params);
        save('NoBc', 'BenchResult');
    case 'OPTIM_BASE_NO_BC'
        Params = OVERWRITE(Params, 'NoBc', 'Base');
        Params.ControlSpec = 'Base';
        Params.TaxSpec = 'Factor';
        Params = COMMON(Params);
        Params.VfiFun = @VFI_FACTOR;
        OPTIM1(Params);
    case 'OPTIM_EXHC_NO_BC'
        Params = OVERWRITE(Params, 'NoBc', 'Base');
        Params.ControlSpec = 'Exhc';
        Params.TaxSpec = 'Factor';
        Params = COMMON(Params);
        Params.VfiFun = @VFI_EXHC_FACTOR;
        OPTIM1(Params);
    case 'CALIBRATE_NO_TYPE'
        Params = OVERWRITE(Params, 'NoType');
        CALIBRATE(Params);
    case 'NO_TYPE'
        Params = OVERWRITE(Params, 'NoType');
        BenchResult = PARTIAL(Params);
        save('NoType','BenchResult');
    case 'OPTIM_BASE_NO_TYPE'
        Params = OVERWRITE(Params, 'NoType', 'Base');
        Params.ControlSpec = 'Base';
        Params.TaxSpec = 'Factor';
        Params = COMMON(Params);
        Params.VfiFun = @VFI_FACTOR;
        OPTIM1(Params);
    case 'OPTIM_EXHC_NO_TYPE'
        Params = OVERWRITE(Params, 'NoType', 'Exhc');
        Params.ControlSpec = 'Exhc';
        Params.TaxSpec = 'Factor';
        Params = COMMON(Params);
        Params.VfiFun = @VFI_EXHC_FACTOR;
        OPTIM1(Params);
    case 'CALIBRATE_NO_IDIO'
        Params = OVERWRITE(Params, 'NoIdio');
        CALIBRATE(Params);
    case 'NO_IDIO'
        Params = OVERWRITE(Params, 'NoIdio');
        BenchResult = PARTIAL(Params);
        save('NoIdio','BenchResult');
    case 'OPTIM_BASE_NO_IDIO'
        Params = OVERWRITE(Params, 'NoIdio', 'Base');
        Params.ControlSpec = 'Base';
        Params.TaxSpec = 'Factor';
        Params = COMMON(Params);
        Params.VfiFun = @VFI_FACTOR;
        OPTIM1(Params);
    case 'OPTIM_EXHC_NO_IDIO'
        Params = OVERWRITE(Params, 'NoIdio', 'Exhc');
        Params.ControlSpec = 'Exhc';
        Params.TaxSpec = 'Factor';
        Params = COMMON(Params);
        Params.VfiFun = @VFI_EXHC_FACTOR;
        OPTIM1(Params);
end

