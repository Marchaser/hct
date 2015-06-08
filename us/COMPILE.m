function COMPILE
% Taxing Human Capital
% @Author: Wenlan Luo
% Compile dpopt with utility functions

make_flag.Define = '-DNON_SEPARABLE';

% make_flag.debug = 'on';
% make_flag.openmp = 'off';
% make_flag.snopt = 'on';
if ~isequal(exist('dpopt_utility','file'),3)
make_flag.donlp2 = 'on';
make_flag.MaxData = 30;
make_flag.MaxDim = 10;
make_dpopt('utility', make_flag);
end

if ~isequal(exist('dpopt_utility_exhc','file'),3)
make_flag.donlp2 = 'on';
make_flag.MaxData = 30;
make_flag.MaxDim = 10;
make_dpopt('utility_exhc', make_flag);
end

if ~isequal(exist('dpopt_utility_factor','file'),3)
make_flag.donlp2 = 'on';
make_flag.MaxData = 30;
make_flag.MaxDim = 10;
make_dpopt('utility_factor', make_flag);
end

if ~isequal(exist('dpopt_utility_factor_trans','file'),3)
make_flag.donlp2 = 'on';
make_flag.MaxData = 30;
make_flag.MaxDim = 10;
make_dpopt('utility_factor_trans', make_flag);
end

if ~isequal(exist('dpopt_utility_exhc_factor','file'),3)
make_flag.donlp2 = 'on';
make_flag.MaxData = 30;
make_flag.MaxDim = 10;
make_dpopt('utility_exhc_factor', make_flag);
end

if ~isequal(exist('dpopt_utility_old','file'),3)
make_flag.snopt = 'off';
make_flag.donlp2 = 'off';
make_flag.MaxData = 30;
make_dpopt('utility_old', make_flag);
end
