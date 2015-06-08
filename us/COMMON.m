function NewParams = COMMON(Params)
% Taxing Human Capital
% @Author: Wenlan Luo
% Prepare common space, random numbers, etc.

% unpack structures
v2struct(Params);

% For replicability
rng(Seed);

% generate A grids
[ATrans, AGrid] = markovappr(0,LogASigma,AGridDispersion,APts);
AGrid = exp(AGrid + LogAMu)';

% epsilon shock
[EpsilonTrans, LogEpsilonGrid] = markovappr(0, LogEpsilonSigma, EpsilonGridDispersion, EpsilonPts);
EpsilonP = EpsilonTrans(1,:);
EpsilonGrid = exp(LogEpsilonGrid + LogEpsilonMu)';
EpsilonCumTrans = cumsum(EpsilonTrans(1,:), 2);
EpsilonCutPoints = floor(NumOfAgents * EpsilonCumTrans);
EpsilonCutPoints(end) = NumOfAgents;  % avoid rounding error
% Note(wenlan): The idea is to generate number of points at each grid point
% consistent with LLN, and then do a random permutation
EpsilonCutPoints = [1 EpsilonCutPoints];
EpsilonShockIdxBase = zeros(1, NumOfAgents);
for i=1:EpsilonPts
    EpsilonShockIdxBase(EpsilonCutPoints(i):EpsilonCutPoints(i+1)) = i;
end
EpsilonShockIdx = zeros(NumOfAgents, Jw);
EpsilonShock = zeros(NumOfAgents, Jw);
for j=1:Jw
    RandIdx = randperm(NumOfAgents);
    EpsilonShockIdx(:,j) = EpsilonShockIdxBase(RandIdx);
    EpsilonShock(:,j) = EpsilonGrid(EpsilonShockIdx(:,j));
end

% generate mesh
[AMesh, KMesh, HMesh] = ndgrid(AGrid, KGrid, HGrid);
ProblemSize = size(KMesh);
NumOfProblems = numel(KMesh);

K1 = zeros(NumOfAgents, 1);

% Note(wenlan): Similarly, generate AIdx exactly satisfying LLN
ACumTrans = cumsum(ATrans(1,:), 2);
ACutPoints = floor(NumOfAgents * ACumTrans);
ACutPoints(end) = NumOfAgents;  % avoid rounding error
ACutPoints = [1 ACutPoints];
% Note(wenlan): I generate h1 based on the conditional distribution h1|a
% for each agent. That is, I first generate a; then compute the conditional
% mean and variance of h1; and finally generate h1 based on the conditional
% mean and variance
ConditionalLogH1Mu = LogH1Mu + LogH1Sigma / LogASigma * HARho * (log(AGrid) - LogAMu);
ConditionalLogH1Variance = (1-HARho^2) * LogH1Sigma^2;
% generate joint initial of a and h1
AIdx = zeros(NumOfAgents, 1);
H1 = zeros(NumOfAgents, 1);
for i=1:1:APts
    AIdx(ACutPoints(i):ACutPoints(i+1)) = i;
    H1(ACutPoints(i):ACutPoints(i+1)) = ConditionalLogH1Mu(i) + randn(ACutPoints(i+1)-ACutPoints(i)+1, 1) * ConditionalLogH1Variance^0.5;
end
AShock = AGrid(AIdx);
H1 = exp(H1);
% random permute of (a, h1)
RandIdx = randperm(NumOfAgents);
AIdx = AIdx(RandIdx);
AShock = AShock(RandIdx);
H1 = H1(RandIdx);

% Note(wenlan): Load benchmark result if it exists in the folder
if ~isequal(exist('BenchResult','var'),1)
    HasBenchResult = exist([BenchResultFileName '.mat'],'file');
    if isequal(HasBenchResult, 2)
        BenchResult = load([BenchResultFileName '.mat']);
        BenchResult = BenchResult.BenchResult;
    end
end
if isequal(exist('BenchResult','var'),1)
    % Process Bench Result
    KBarIdx = 1:KShrinkFactor:length(KGrid);
    HBarIdx = 1:HShrinkFactor:length(HGrid);
    KBarPts = length(KBarIdx);
    HBarPts = length(HBarIdx);
    KBarGrid = KGrid(KBarIdx);
    HBarGrid = HGrid(HBarIdx);
    % Bench Policy has the memory order (age, alpha, KBar, HBar)
    BenchSPolicy = permute(squeeze(BenchResult.VfiResult.PolicyYoung(1,:,KBarIdx,HBarIdx,:)), [4 1 2 3]);
    BenchNPolicy = permute(squeeze(BenchResult.VfiResult.PolicyYoung(2,:,KBarIdx,HBarIdx,:)), [4 1 2 3]);
    BenchKpPolicy = permute(squeeze(BenchResult.VfiResult.PolicyYoung(3,:,KBarIdx,HBarIdx,:)), [4 1 2 3]);
    BenchSimulationS = BenchResult.SimulateResult.s;
    BenchSimulationN = BenchResult.SimulateResult.n;
    BenchSimulationK = BenchResult.SimulateResult.k;
    BenchSimulationH = BenchResult.SimulateResult.h;
end

r = Gamma * KLRatio^(Gamma-1) - Delta;
w = (1-Gamma) * KLRatio^Gamma;

% pack structure
clear Params;
NewParams = v2struct();
end
