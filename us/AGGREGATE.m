function [AggregateResult, ExitFlag] = AGGREGATE(SimulateResult, Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Compute aggregates based on simulated distribution

% unpack structure
v2struct(Params);
v2struct(SimulateResult);

% if need to print information

if (isfield(Params, 'AggregateDisplay'))
    AggregateDisplay = Params.AggregateDisplay;
else
    AggregateDisplay = true;
end

% compute per capita aggregates
C = WMean(c(:), population(:));
TaxC = C * TauC;
K = WMean(k(:), population(:));
L = WMean(h(:) .* n(:), population(:));
% Tss distributed to retirement
Tss = sum(sum(tss(:,1:Jw).*population(:,1:Jw))) / ...
    sum(sum(population(:,Jw+1:end)));
TaxInc = WMean(tax(:), population(:));
Tax = TaxInc + TaxC;
G = Tax;
% accidental bequest is the saving by dead agents
Tr = sum(sum(k(:,2:J).*(1-survival(:,1:J-1)).*population(:,1:J-1))) / ...
    sum(population(:));
Y = K^Gamma * L^(1-Gamma);
GYRatio = G / Y;
KYRatio = K / Y;
KLRatio = K / L;
r = Gamma * KLRatio^(Gamma-1) - Delta;
w = (1-Gamma) * KLRatio^Gamma;
MeanHours = WMean(reshape(n(:,1:Jw)+s(:,1:Jw),1,[]),reshape(population(:,1:Jw),1,[]));
% MeanHours = WMean(reshape(n(:,1:Jw),1,[]),reshape(population(:,1:Jw),1,[]));
% wage profile
VarLogWage = std(log(h(:,1:40))).^2;
MeanWage = mean(h(:,1:40));
MeanMedianRatio = mean(h(:,1:40)) ./ median(h(:,1:40));
WageGini = ginicoeff(h(:,1:40));
% welfare
V = mean(v);
CumulativeSurvivalRate = cumprod([1 SurvivalRate(1:J-1)']);
CumulativeBeta = Beta .^ [0:J-1];
U = repmat(CumulativeBeta .* CumulativeSurvivalRate, NumOfAgents, 1) .* ...
    ((c.^Chi.*l.^(1-Chi)).^(1-Sigma1)) / (1-Sigma1);
Welfare = mean(sum(U, 2));
% Welfare = V;

MomentsModel = [MeanWage/MeanWage(1) VarLogWage MeanMedianRatio WageGini];

Stats = [];
if (AggregateDisplay == true)
% sumary statistics
% working age
variableList = {'k','h','s','n','c','l','tax','taxk','taxn'};
StatsList = {'max', 'min', 'mean', 'std'};
fprintf('%-8s', '1:Jw');
for j=1:length(StatsList)
    fprintf('%8s', StatsList{j});
end
fprintf('\n');
for i=1:length(variableList)
    if (isfield(SimulateResult, variableList{i}))
        fprintf('%-8s', variableList{i});
        for j=1:length(StatsList)
            eval(['Stats.' variableList{i} StatsList{j} '=' StatsList{j} '(reshape(' variableList{i} '(:,1:Jw),1,[]));']);
            eval(['fprintf(''%8.4f'',Stats.' variableList{i} StatsList{j} ')']);
        end
        fprintf('\n');
    end
end

variableList = {'k','c','tax','ktax','ntax'};
StatsList = {'max', 'min', 'mean', 'std'};
fprintf('%-8s', '1:J');
for j=1:length(StatsList)
    fprintf('%8s', StatsList{j});
end
fprintf('\n');
for i=1:length(variableList)
    if (isfield(SimulateResult, variableList{i}))
        fprintf('%-8s', variableList{i});
        for j=1:length(StatsList)
            eval(['Stats.' variableList{i} StatsList{j} '=' StatsList{j} '(reshape(' variableList{i} '(:,1:J),1,[]));']);
            eval(['fprintf(''%8.4f'',Stats.' variableList{i} StatsList{j} ')']);
        end
        
        fprintf('\n');
    end
end

fprintf('K, L, Y, G, G/Y, K/Y:, %g, %g, %g, %g, %g, %g\n', K, L, Y, G, GYRatio, KYRatio);
fprintf('r, w, Tss, Tr, Tax_inc, mean(hours):, %g, %g, %g, %g, %g, %g\n', r, w, Tss, Tr, TaxInc, MeanHours);
fprintf('var(log_wage), mean(wage): %g, %g, %g\n', VarLogWage(1), VarLogWage(end), max(MeanWage)/MeanWage(1));
fprintf('V, Welfare: %g, %g\n', V, Welfare);

end

% return
AggregateResult = v2struct(C, TaxC, K, L, Tss, TaxInc, Tax, G, Tr, Y, GYRatio, KYRatio, KLRatio, r, w, MeanHours, VarLogWage, MeanWage, MeanMedianRatio, WageGini, V, Welfare, MomentsModel, Stats);
ExitFlag = 0;
end

function s = WMean(x, weight)
s = sum(x.*weight)./sum(weight);
end