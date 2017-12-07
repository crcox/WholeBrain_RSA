subjix = 1:23;
cvholdout = 0;
PERMUTATION_INDEX{1} = ones(100,2);
BIAS = true;
LambdaSeq = 'Inf';
normalize = 'zscore';
regularization = 'GROWL2';
lambda = rand(1,24);
lambda1 = rand(1,24);

AdlasInstances = AdlasContainer( ...
    'subject', subjix, ...
    'RandomSeed', 1:size(PERMUTATION_INDEX{1},2), ...
    'cvholdout', cvholdout, ...
    'bias', BIAS, ...
    'LambdaSeq', LambdaSeq, ...
    'normalize', normalize, ...
    'regularization', regularization, ...
    'HYPERBAND', struct('lambda', lambda, 'lambda1', lambda1));

n = numel(AdlasInstances);
for i = 1:n
    AdlasInstances(i).Adlas.testError = randn(1);
end

[~,ix] = hyperband_pick_top_n(AdlasInstances, 4);
disp(AdlasInstances(ix));
