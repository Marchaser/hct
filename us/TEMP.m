
max(reshape(abs(VfiTransResult.PolicyOld(:,:,1)-VfiTransResult.PolicyOld(:,:,2)),1,[]))
max(reshape(abs(VfiTransResult.PolicyOld(:,:,3)-VfiTransResult.PolicyOld(:,:,2)),1,[]))

VfiTransResult.PolicyOld(:,:,1)
VfiTransResult.PolicyOld(:,:,2)

VfiTransResult.ValueOld(:,:,6)
max(reshape(abs(VfiTransResult.ValueOld(:,:,6) - VfiTransResult.ValueOld(:,:,5)),1,[]))
VfiTransResult.ValueOld(:,:,6) - VfiTransResult.ValueOld(:,:,5)

max(reshape(abs(VfiTransResult.PolicyYoung(:,:,:,:,:,5)-VfiTransResult.PolicyYoung(:,:,:,:,:,4)),1,[]))

max(reshape(abs(VfiTransResult.ValueYoung(:,:,:,:,6)-VfiTransResult.ValueYoung(:,:,:,:,5)),1,[]))
VfiTransResult.ValueYoung(:,:,:,:,6)-VfiTransResult.ValueYoung(:,:,:,:,5)

load('SimulateTransFullResultBase');
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

max(reshape(abs(VfiSs1Result.ValueYoung-VfiTransResult.ValueYoung(:,:,:,:,:,1)),1,[]))

