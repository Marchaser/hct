function NewParams = MOMENTS(Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Attach data moments to Params
v2struct(Params);

GYRatioTarget = 0.17;
KYRatioTarget = 2.7;
MeanHoursTarget = 0.33;

PsidStats = csvread('psid_stats.csv',1,0);

PsidStats = PsidStats(1:40, :);

MeanWagePsid = PsidStats(:, 1);
VarLogWagePsid = PsidStats(:, 2);

MeanWagePsid = MeanWagePsid / MeanWagePsid(1);

MeanMedianRatioPsid = PsidStats(:,3);
WageGiniPsid = PsidStats(:,4);

MomentsData = [ ...
    MeanWagePsid' ...
    VarLogWagePsid' ...
    MeanMedianRatioPsid' ...
    WageGiniPsid' ...
    ];

MomentsWeight = [ones(1,40) 5*ones(1,40)];
% MomentsWeight = ones(size(MomentsData));

clear Params;
NewParams = v2struct;
end