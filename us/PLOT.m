function PLOT(AggregateResult, Params)
ModelVsDataFigure = figure;
subplot(2,2,1);
hold on;
grid on;
plot([21:60],AggregateResult.MeanWage / AggregateResult.MeanWage(1),'k--','LineWidth',1);
plot([21:60],Params.MeanWagePsid,'b-','LineWidth',1);
xlabel('Age');
title('Mean of wage');
legend('Model', 'Data', 'Location','SouthEast');
hold off;

subplot(2,2,2);
hold on;
grid on;
plot([21:60],AggregateResult.VarLogWage,'k--','LineWidth',1);
plot([21:60],Params.VarLogWagePsid,'b-','LineWidth',1);
xlabel('Age');
title('Variance of log wage');
legend('Model', 'Data', 'Location','SouthEast');
hold off;

subplot(2,2,3);
hold on;
grid on;
plot([21:60],AggregateResult.WageGini,'k--','LineWidth',1);
plot([21:60],Params.WageGiniPsid,'b-','LineWidth',1);
xlabel('Age');
title('Gini of wage');
legend('Model', 'Data', 'Location','SouthEast');
hold off;

subplot(2,2,4);
hold on;
grid on;
plot([21:60],AggregateResult.MeanMedianRatio,'k--','LineWidth',1);
plot([21:60],Params.MeanMedianRatioPsid,'b-','LineWidth',1);
xlabel('Age');
title('Skewness of wage (mean/median)');
legend('Model', 'Data', 'Location','SouthEast');
hold off;

set(findall(gcf,'type','text'),'FontSize',14,'FontName','Times New Roman')
print(ModelVsDataFigure, 'graph/ModelVsData.pdf', '-dpdf');

end