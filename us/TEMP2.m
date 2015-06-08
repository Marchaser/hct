%{
for j=1:Params.Jw
    for iHGroup = 1:10 
        idx = SimulateSs0FullResult.h(:,j)<quantile(SimulateSs0FullResult.h(:,j),0.1*(iHGroup)) & SimulateSs0FullResult.h(:,j)>=quantile(SimulateSs0FullResult.h(:,j),0.1*(iHGroup-1));
        SupportByJH(iHGroup,j) = mean(vDiff(idx,j)>0 ...
            .* SimulateSs0FullResult.population(idx,j));
    end
end

[HGroupMesh,JwMesh] = ndgrid([1:10],[1:Params.Jw]);
mesh(HGroupMesh,JwMesh,SupportByJH);
%}

SupportByJK = [];
for j=1:Params.J-2
    for iKGroup = 1:10 
        idx = SimulateSs0FullResult.k(:,j)<=quantile(SimulateSs0FullResult.k(:,j),0.1*(iKGroup)) & SimulateSs0FullResult.k(:,j)>=quantile(SimulateSs0FullResult.k(:,j),0.1*(iKGroup-1));
        SupportByJK(iKGroup,j) = mean(vDiff(idx,j)>0 ...
            .* SimulateSs0FullResult.population(idx,j));
    end
end

[KGroupMesh,JwMesh] = ndgrid([0.1:0.1:1],[1:Params.J-2]);
figure;
surf(KGroupMesh,JwMesh,SupportByJK);
xlabel('Capital stock percentiles');
ylabel('Age');
title('Political Support of transition');
shading flat;
colormap(gray);
view(90,90);
% print(gcf, 'graph/support_base.pdf', '-dpdf');
print(gcf, 'graph/support_exhc.pdf', '-dpdf');
