function ANALYZE

format short

Bench = load('BenchResult');
Bench = Bench.BenchResult;
Base = load('BaseAtOptim');
Base = Base.BaseAtOptim;
BaseSub = load('BaseAtSubOptim');
BaseSub = BaseSub.BenchResult;

ExHc = load('ExhcAtOptim');
ExHc = ExHc.ExhcAtOptim;

%{
ExN = load('ExnAtOptim');
ExN = ExN.ExnAtOptim;
%}

inspect_tax_function(Bench, Base);
% inspect_tax_function(Bench, ExHc);
% compute_tax_stats(Bench,Bench);
% compute_tax_stats(Base,Bench);
% compute_tax_stats(ExHc,Bench);

%{
compute_tax_stats(Bench);
compute_tax_stats(Base);
compute_tax_stats(ExHc);
compute_tax_stats(ExN);

inspect_borrowing_constraints(Bench);
inspect_borrowing_constraints(Base);
inspect_borrowing_constraints(ExHc);
%}


Bench = append_structures(Bench);
% bench_vs_optimal_plot(Bench,Base);
Base = append_structures(Base); 
% bench_vs_optimal_plot(Bench,Base);

bench_vs_optimal(Bench,Base);
bench_vs_optimal(Bench,Base);
%{
ExHc = append_structures(ExHc);
bench_vs_optimal(Bench,ExHc);
%}
BaseSub = append_structures(BaseSub);
bench_vs_optimal(Bench,BaseSub);

% model_performance(Bench);
end

function bench_vs_optimal_plot(Bench,Optimal)
Params = COMMON(Bench.Params);

Y = Optimal.Y;
current_Y = 51;
% adjust model scale to real scale
fileList = {'Bench', 'Optimal'};
variableList = {'c', 'k', 'tax'};
for i=1:length(variableList)
    variable = variableList{i};
    for j=1:length(fileList)
        file = fileList{j};
        eval(sprintf('%s.%s=%s.%s / Y * current_Y;', file, variable, file, variable));
    end
end

fileList = {'Bench', 'Optimal'};
fileDescList = {'Bench', 'Optimal'};
% variableList = {'c', 'k', 'n', 's','h','tax'};
variableList = {'c', 'k', 'n', 's'};
titleList = {'Consumption', 'Assets', 'Work hours', 'Human capital investment', 'Human capital stock', 'Income taxes'};
ylabelList = {'Consumption ($1,000)', 'Assets ($1,000)', 'Work hours', 'Human capital investment', 'Human capital stock', 'Income taxes ($1,000)'};
yrange = {[0 70], [0 500], [0 0.35], [0 0.16], [0 2], [0 18]};
xrange = [20 100];
% variableList{5:6}={};
SubPlotLocation = {[1 1], [1 2], [2 1], [2 2]};
OptimalVsBenchFigure = figure;
for i=1:length(variableList)
    variable = variableList{i};
    plotList = [];
    for j=1:length(fileList)
        file = fileList{j};
        eval(sprintf('mean_%s%d_low = mean(%s.%s(Params.AIdx==1,:));', variable, j, file, variable));
        eval(sprintf('mean_%s%d_high = mean(%s.%s(Params.AIdx==2,:));', variable, j, file, variable));
        eval(sprintf('plotList = [plotList;mean_%s%d_low;mean_%s%d_high];', variable, j, variable, j));
    end
    % eval(sprintf('%sFigure = figure', variable));
    
    subplot(2,2,i);
    hold on;
    grid on;
    plot([1:size(plotList,2)]+20, plotList(1,:), '-k', 'LineWidth', 2);
    plot([1:size(plotList,2)]+20, plotList(3,:), '--k', 'LineWidth', 2);
    plot([1:size(plotList,2)]+20, plotList(2,:), '-b', 'LineWidth', 2);
    plot([1:size(plotList,2)]+20, plotList(4,:), '--b', 'LineWidth', 2);
    legend(fileDescList);
    xlabel('Age');
    ylim(yrange{i});
    xlim(xrange);
    eval(sprintf('ylabel(''%s'');', ylabelList{i}));
    eval(sprintf('title(''%s'');', titleList{i}));
%     eval(sprintf('set(findall(%sFigure,''type'',''text''),''FontSize'',24,''FontName'',''Times New Roman'');', variable));
    set(gca,'fontsize',18);
    set(gca,'fontname','Times New Roman');
    % eval(sprintf('print(%sFigure, ''graph/by_a_%s.pdf'', ''-dpdf'');', variable, variable));
    hold off;
    % eval(sprintf('close(%sFigure);', variable));
end
set(findall(OptimalVsBenchFigure,'type','text'),'FontSize',12,'FontName','Times New Roman');
print(OptimalVsBenchFigure, 'graph/OptimalVsBench', '-dpdf');
end

function bench_vs_optimal(Bench,Optimal)
v2struct(Bench.Params);
DeltaH = (sum(sum(Optimal.h(:,1:Jw) .* Optimal.population(:,1:Jw))) / sum(sum(Bench.h(:,1:Jw) .* Bench.population(:,1:Jw)))- 1) *100;
DeltaL = (Optimal.L / Bench.L - 1) * 100;
DeltaK = (Optimal.K / Bench.K - 1) * 100;
DeltaY = (Optimal.Y / Bench.Y - 1) * 100;

DeltaN = (sum(sum(Optimal.n(:,1:Jw) .* Optimal.population(:,1:Jw))) / sum(sum(Bench.n(:,1:Jw) .* Bench.population(:,1:Jw))) - 1) *100;
DeltaLeisure = (sum(sum(Optimal.l(:,1:J) .* Optimal.population(:,1:J))) / sum(sum(Bench.l(:,1:J) .* Bench.population(:,1:J))) - 1) *100;
DeltaS = (sum(sum(Optimal.s(:,1:Jw) .* Optimal.population(:,1:Jw))) / sum(sum(Bench.s(:,1:Jw) .* Bench.population(:,1:Jw))) - 1) *100;
DeltaC = (sum(sum(Optimal.c(:,1:J) .* Optimal.population(:,1:J))) / sum(sum(Bench.c(:,1:J) .* Bench.population(:,1:J))) - 1) *100;

DeltaLByN = (sum(sum(Optimal.n(:,1:Jw) .* Bench.h(:,1:Jw) .* Optimal.population(:,1:Jw))) / sum(sum(Bench.n(:,1:Jw) .* Bench.h(:,1:Jw) .* Bench.population(:,1:Jw))) - 1) *100;
DeltaLByH = (sum(sum(Bench.n(:,1:Jw) .* Optimal.h(:,1:Jw) .* Optimal.population(:,1:Jw))) / sum(sum(Bench.n(:,1:Jw) .* Bench.h(:,1:Jw) .* Bench.population(:,1:Jw))) - 1) *100;

Cev = ((Optimal.Welfare / Bench.Welfare)^(1/(Chi*(1-Sigma1))) - 1)*100;

display([DeltaL DeltaLByN DeltaLByH]);
display([DeltaLByN/DeltaL DeltaLByH/DeltaL 1-(DeltaLByN+DeltaLByH)/DeltaL]);
display([DeltaN DeltaS DeltaH DeltaL DeltaK DeltaY DeltaC DeltaLeisure Cev]');
end

function NewEq = append_structures(Eq)
v2struct(Eq.SimulateResult);
v2struct(Eq.AggregateResult);
v2struct(Eq.Params);
Params = Eq.Params;
clear Eq;
NewEq = v2struct;
end

function model_performance(Eq)
v2struct(Eq.Params);
v2struct(Eq.SimulateResult);

% inequality measure
Wealth = reshape(k(:,1:J), 1, []);
WealthGini = ginicoeff(Wealth)
Earnings = reshape(ninc(:,1:Jw), 1, []);
EarningsGini = ginicoeff(Earnings)
Income = kinc(:) + ninc(:);
IncomeGini = ginicoeff(Income)

% compare with earnings quantile share
% from Diaz-Gimenez et al. (2011),
EarningsShareData = [-0.1 4.2 11.7 20.8 63.5];
Perc = [0 0.2 0.4 0.6 0.8 1.0];
EarningsQuantile = quantile(Earnings, Perc);
EarningsTotal = [];
for i=1:length(Perc)-1
    EarningsTotal(i) = sum(((Earnings < EarningsQuantile(i+1) & (Earnings>EarningsQuantile(i))).*Earnings));
end
EarningsShareModel = EarningsTotal / sum(Earnings) * 100;
display(EarningsShareData);
display(EarningsShareModel);

WealthShareData = [-0.2 1.1 4.5 11.2 83.4];
Perc = [0 0.2 0.4 0.6 0.8 1.0];
WealthQuantile = quantile(Wealth, Perc);
WealthTotal = [];
for i=1:length(Perc)-1
    WealthTotal(i) = sum(((Wealth < WealthQuantile(i+1) & (Wealth>WealthQuantile(i))).*Wealth));
end
WealthShareModel = WealthTotal / sum(Wealth) * 100;
display(WealthShareData);
display(WealthShareModel);

IncomeShareData = [2.8 6.7 11.3 18.3 60.9];
Perc = [0 0.2 0.4 0.6 0.8 1.0];
IncomeQuantile = quantile(Income, Perc);
IncomeTotal = [];
for i=1:length(Perc)-1
    IncomeTotal(i) = sum(((Income < IncomeQuantile(i+1) & (Income>IncomeQuantile(i))).*Income));
end
IncomeShareModel = IncomeTotal / sum(Income) * 100;
display(IncomeShareData);
display(IncomeShareModel);
end

function inspect_tax_function(Bench, Optimal)
v2struct(Optimal.Params);

% Shape of GS function
% Get income from bench
Y = Bench.AggregateResult.Y;
Y2004 = 42;
% income from 0:2Y
Inc = 0:0.1:2*Y2004;
% scale to model
IncModel = Inc / Y2004 * Y;
% apply to tax function
TaxModel = GS(IncModel, Tau0N, Tau1N, Tau2N);
TauMModel = GS_PRIME(IncModel, Tau0N, Tau1N, Tau2N);
TauModel = TaxModel ./ IncModel;
% plot(IncModel, TauModel);
figure;
plot(IncModel, TauMModel);

% find deductibles
HalfMarginalRate = Tau0N/2;
% solve income that gives half marginal rate
IncAtHalfMarginalRate = fzero(@(x) GS_PRIME(x,Tau0N,Tau1N,Tau2N)-HalfMarginalRate, 0.1);
% convert income to norminal
Deductibles = IncAtHalfMarginalRate / Y * Y2004
end

function inspect_borrowing_constraints(Eq)
v2struct(Eq.Params)
v2struct(Eq.SimulateResult);
display(sum(k==KMin));
end

function compute_tax_stats(Eq,Bench)
v2struct(Eq.Params);
v2struct(Bench.SimulateResult);
inc = ninc(:,1:Jw) + kinc(:,1:Jw);
inc = inc(:);

switch TaxSpec
    case 'Total'
        compute_stats_given_inc_tax_and_tau(inc(:), Tau0Gs, Tau1Gs, Tau2Gs);
    case 'Factor'
        display('Capital income tax:');
        compute_stats_given_inc_tax_and_tau(kinc(:), Tau0K, Tau1K, Tau2K);
        display('Labor income tax:');
        compute_stats_given_inc_tax_and_tau(reshape(ninc(:,1:Jw),[],1), Tau0N, Tau1N, Tau2N);
end
end

function compute_stats_given_inc_tax_and_tau(Inc, Tau0, Tau1, Tau2)
MedianInc = median(Inc);
MedianIncHalf = MedianInc*0.5;
MedianIncByTwo = MedianInc*2;

marginal_tax_rate = @(x) GS_PRIME(x,Tau0,Tau1,Tau2);

TaxRateAtMedianInc = marginal_tax_rate(MedianInc);
TaxRateAtMedianIncHalf = marginal_tax_rate(MedianIncHalf);
TaxRateAtMedianIncByTwo = marginal_tax_rate(MedianIncByTwo);

PwHalfTwo = 1 - (1-TaxRateAtMedianIncByTwo)/(1-TaxRateAtMedianIncHalf);
PwHalfOne = 1 - (1-TaxRateAtMedianInc)/(1-TaxRateAtMedianIncHalf);

Atr = mean(marginal_tax_rate(Inc));

display(Atr);
display([PwHalfOne PwHalfTwo]);
end
