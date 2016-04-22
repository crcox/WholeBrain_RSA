% generate data and similarity structure

% First, let's take a very simple example with a 25 voxel brain. There are
% several important concepts at play. Let's imagine a study where participants
% made semantic judgements about 50 items. We have an estimate of the semantic
% similarity structure derived from, e.g., a set of feature norms. Each item is
% associated with, say, 300 features that can be either present or absent. The
% similarity among these feature norms is summarized by their pairwise euclidean
% distances, which will yield a 50x50 symmetric matrix with zeros on the diagonal.

% That is the scenario you should have in mind, but for the purposes of the
% simulated example we will actually generate the data so there is a known
% relationship between the ``MRI data'' (X) and the similarity structure (S) by
% generating the S from X via a set of weights, U.

% Rather than being a vector with lenght nVoxels, U is a nVoxel x nDimension
% matrix. Dimension here is related to the ``rank'' of S, where structured
% solutions will have lower rank. Put another way, in a singular value
% decomposition (e.g. principle components analysis) of a matrix, the matrix
% will be decomposed into as many components as their are features, but only a
% small number of them will have large eigenvalues. When you plot the
% eigen-values of a non-random matrix, you will typically see that the
% magnitude of the eigen-values drops precipitously after the first few. This
% means that most of the variance in the matrix can be explained in a few
% dimensions, which is another way of saying there is ``similarity structure''
% expressed within the matrix.

nItems = 50;
nVoxels = 25;
nDimensions = 3;
nSignificantVoxels = 25;

% U will encode the relationship between X and sqrt(S)=Y. This is a little
% trick. Y=sqrt(S), and U=sqrt(W). The relationship between S and X would be
% modeled as S=XWX', but this is a difficult problem to solve. Y=XU is much
% easier to solve, and is related via a sqrt to the solution we want and so
% this is how we do it in practice.
X = randn(nItems,nVoxels);
U = zeros(nVoxels, nDimensions);
U(1:nSignificantVoxels,:) = randn(nSignificantVoxels, nDimensions);
W = U*U';

% N.B. Y will have as many columns as U, and the number of columns corresponds directly to the ``rank'' of S, because S=YY'
Y = X * U;
S = Y*Y';

% Now we have a ground truth for S, and W, and we have X. In practice, we will
% have only S and X and we want to solve for W. Because we have more items than
% voxels in this toy example, we can solve this relatively simply.

% The following closed form solution for the weight estimates, Uz, is based on
% the matrix formalization of regression. Note that the `z' suffix will mark
% estimated values. The `1' suffix means that the predictions or values are
% about the un-trained holdout set. We will fit these parameters, holding out
% 20% of the items for a test set.
c = cvpartition(50, 'HoldOut', .2)
X1 = X(test(c),:);
Y1 = Y(test(c),:);
X2 = X(training(c),:);
Y2 = Y(training(c),:);
Uz = pinv(X1) * Y1;

% The following is one possible way to compute an error metric. It's important
% to note that this error metric, in this context, does not have a natural
% value associated with ``chance'' or zero-structure, and it should be
% determined (in practice) by permutation test. Note also that, in practice, we
% should have done cross validation.
S1 = Y1 * Y1';
Yz1 = X(test(c),:) * Uz;
Sz1 = Yz1 * Yz1';
err_Y = norm(Y1 - Yz1,'fro')/norm(Y1,'fro');
err = norm(S1 - Sz1,'fro')/norm(S1,'fro');

% err_Y and err are fundamentally related, but both are computed for
% completeness. Error is essentially zero.

% You might be wondering how this differs from a standard RSA, where the
% similarity structure in the MRI data is no modeled, but assessed directly. We
% can make this comparison explicit by simply computing the euclidean distance
% among the rows of X, and seeing who closely this structure corresponds to the
% true S. Here, I am computing a ``similarityStructureScore'', which is simply
% the sum over the element-wise products between the two matrices. Larger
% scores are better---as the similarity between the two matrices increases, so
% does the score. Imporant to note, again, that this is on a totally arbitrary
% scale, and so can only be interpretted relative to something else, or to
% permutation tests.
d1 = pdist(X1);
D1 = squareform(d1);
similarityStructureScore = zeros(1,2);
similarityStructureScore(1) = S1(:)' * D1(:);
similarityStructureScore(2) = S1(:)' * Sz1(:);

% It is clear that the similarity score between S and Sz is much, *much* higher
% than S and D. What is going on? In a standard RSA, there is no model. It
% assumes that the similarity structure is expressed equally over all voxels
% considered, and that all voxels are equally relevant to all underlying
% dimensions of the similarity space. Network RSA allows the relevance of each
% voxel to each dimension to be weighted, and the predicted similarity
% structure is related to the weighted sum of neural activity. By anology,
% consider multivoxel (pattern) classification analysis. In MVPA, the voxel
% activity is modeled to make category predictions, and does not assume that
% all voxels are equally important. Network RSA brings this flexibility to RSA.

% Another way of driving this point home is to realize that the standard RSA
% analysis is similar to assuming that the weight on all voxels is 1. With this
% in mind, we can "make predictions" using the "standard model", which is
% simply the identity matrix.
Si1 = X1*eye(nVoxels)*X1';
err = zeros(1,2);
err(1) = norm(S1 - Si1,'fro')/norm(S1,'fro');
err(2) = norm(S1 - Sz1,'fro')/norm(S1,'fro');

% This reiterates the comparison above: Standard RSA suggests there is no
% structure in the data, but Network RSA shows that the data actually express
% the structure nearly perfectly.

% At this point, many of the core concepts have been laid out and compared with
% more basic RSA approaches. Next, we consider another key aspect of Network RSA,
% namely that sparse signal can be identified among a large set of features.
% Put another way, this means that Network RSA can discover a subset of voxels
% within a large volume (perhaps all of cortex) that encode similarity
% structure.

nItems = 50;
nVoxels = 1000;
nDimensions = 3;
nSignificantVoxels = 25;

X = randn(nItems,nVoxels);
U = zeros(nVoxels, nDimensions);
U(1:nSignificantVoxels,:) = randn(nSignificantVoxels, nDimensions) * 2;
W = U * U';
Y = X * U;
S = Y * Y';

c = cvpartition(50, 'HoldOut', .2)
X1 = X(test(c),:);
Y1 = Y(test(c),:);
S1 = Y1 * Y1';
X2 = X(training(c),:);
Y2 = Y(training(c),:);
Uz = pinv(X2) * Y2;

Yz1 = X1 * Uz;
Sz1 = Yz1 * Yz1';
err_Y = norm(Y1 - Yz1,'fro')/norm(Y1,'fro');
err = norm(S1 - Sz1,'fro')/norm(S1,'fro');

% Now, the error is very high. How high? We would need to do a permutation test
% to have a point of reference.
nPerm = 1000;
perm_err = zeros(nPerm,1);
for iPerm = 1:nPerm
  Uz = pinv(X2(randperm(c.TrainSize),:)) * Y2;
  Yz1 = X1 * Uz;
  Sz1 = Yz1 * Yz1';
  perm_err(iPerm) = norm(S1 - Sz1,'fro')/norm(S1,'fro');
end
fprintf('Proportion of permutations with error less than CV error: %.3f\n',nnz(perm_err < err)/nPerm)

% Clearly, we are unable to do better than chance when the signal is contained
% in a sparse subset of voxels. The problem is that there is no unique
% solution, and the training set is being over-fit. One solution to the
% overfitting problem is LASSO. In this case, because there are multiple
% columns of W, the solution is actually to use Group LASSO. Group Lasso is
% implemented within the WholeBrain_RSA toolbox.
addpath('../src');

% LASSO involves a free parameter that scales how much the model is penalized
% for having non-zero weights. If you want matlab to find a good value for
% this, run the line before. Otherwise, you can just try a few values and see
% if one works. In production you'll want to do some sort of principled
% parameter search, but for now something in the 0--10 range will probably do.
[lam1, err_L1L2] = fminbnd(@(x) optimizeAdlas1(X,Y,c,x), 0, 10)

[Uz, info] = Adlas1(X2, Y2, lam1);

% This error looks better, but we'll need to do a new permutation test to be sure.
% This is a slower procedure, so I'm going to dial back the permutation count...
nPerm = 10;
perm_err_L1L2 = zeros(nPerm,1);
for iPerm = 1:nPerm
  disp(iPerm)
  [Uz, info] = Adlas1(X2, Y2, lam1);
  Uz = pinv(X2(randperm(c.TrainSize),:)) * Y2;
  Yz1 = X1 * Uz;
  Sz1 = Yz1 * Yz1';
  perm_err_L1L2(iPerm) = norm(S1 - Sz1,'fro')/norm(S1,'fro');
end
fprintf('Proportion of permutations with error less than CV error: %.3f\n',nnz(perm_err_L1L2 < err_L1L2)/nPerm)

% This permutation test suggests that the difference between the predicted and
% true similarity structure is smaller than would be expected by chance.

% One limitation of LASSO (and Group LASSO) is that is seeks the sparsest
% possible solution, and so if there are multiple correlated features, LASSO
% will select one and drop the others. Ordered Weighted LASSO (OWL) is an
% alternative that helps retain correlated voxels and so makes more sense in a
% neuroscience context. Group OWL (GrOWL) is implemented in the WholeBrain_RSA
% package. Of course, these simulated data are generated at random and so we
% shouldn't expect an increase, and perhaps may even see a decrease, in
% preformance on these i.i.d. simulated data.

% GrOWL (implemented in Adlas2) takes 2 lambda arguments. The first affects how
% much to care about correlated voxels, and the second effects how much to care
% about sparsity.
[lamvals, err_growl] = fminsearch(@(x) optimizeAdlas2(X,Y,c,x), [.001,1])

[Uz, info] = Adlas2(X2, Y2, lamvals(1), lamvals(2));

% Probably due to the i.i.d. nature of the data, this yields a high error.
% Again, though, we will not have context for this error term without a
% permutation test.
nPerm = 10;
perm_err_growl = zeros(nPerm,1);
for iPerm = 1:nPerm
  disp(iPerm)
  [Uz, info] = Adlas2(X2, Y2, lam1);
  Uz = pinv(X2(randperm(c.TrainSize),:)) * Y2;
  Yz1 = X1 * Uz;
  Sz1 = Yz1 * Yz1';
  perm_err_growl(iPerm) = norm(S1 - Sz1,'fro')/norm(S1,'fro');
end
fprintf('Proportion of permutations with error less than CV error: %.3f\n',nnz(perm_err_L1L2 < err_L1L2)/nPerm)

% We have now introduced Network RSA, and how it can be achieved with Group
% LASSO and GrOWL using functions in the WholeBrain_RSA package. However, the
% computational intensity of working with these methods may have impressed
% themselves upon you at this point. Parameter searching, k-fold cross
% validation, and permutation testing, mean looping over relatively long
% running processes thousands of times. For this reason, WholeBrain_RSA was
% written with distribued computing in mind.

% The data
% ========
% The fMRI data should be formatted in a time x voxel matrix, X. Each row of
% this matrix is a training example, and so should include all voxels that the
% model may be trained on/evaluated with. By ``time'', I also mean ``item''.
% That is, the method does not require that the data is the raw time series. In
% fact, in my experience it is more common to first peform an item-wise
% deconvolution which will result in a single volume of beta weights for each
% item. In that case, the item x voxel matrix will contain the fitted betas.

% For example, imagine a study with 100 unique items, sampled equally from two
% categories or belonging to two experimental conditions. Further, imagine that
% there are 10,000 voxels in the cortex of this subject.
nItems = 100;
nVoxels = 10000;
X = randn(nItems, nVoxels);

% Let S represent the target similarity structure. This might correspond to
% conditions or items, but the number of rows and columns in S must match the
% number of rows in X. This might involve averaging rows of X together.
Y = struct('visual',randn(100,3),'semantic',randn(100,8));
VISUAL = squareform(pdist(Y.visual));
SEMANTIC = squareform(pdist(Y.semantic));

% The metadata
% ============
% Targets
% -------
% Information about targets (i.e., possible y vectors) should be stored in a
% structure with 5 required fields: 'label', 'similarity', 'sim_source', 'sim_metric', and 'target'.
TARGETS(1).label = 'visual';
TARGETS(1).type = 'similarity'; % {'category','similarity'}
TARGETS(1).sim_source = 'chamfer';
TARGETS(1).sim_metric = 'chamfer';
TARGETS(1).target = VISUAL;

TARGETS(2).label = 'semantic';
TARGETS(2).type = 'similarity'; % {'category','similarity'}
TARGETS(2).sim_source = 'featurenorms';
TARGETS(2).sim_metric = 'cosine';
TARGETS(2).target = SEMANTIC;

% Cross-validation
% ----------------
% The methods in WholeBrain_RSA will not generate CV indexes for you. This is
% to help promote replicability of results. So you will need to provide a cross
% validation scheme ahead of time. If you are concerned about results being
% specific to a particular cross-validation scheme, you can specify multiple
% schemes.
%
% For example, let's set up 10 cross validation schemes.
nschemes = 10;
nfolds = 10;
SCHEMES = zeros(nItems, nschemes);
for iScheme = 1:nschemes
  c = cvpartition(nItems,'KFold', nfolds);
  for iFold = 1:nfolds
    SCHEMES(:,iScheme) = SCHEMES(:,iScheme) + (test(c, iFold) * iFold);
  end
end

% Filters
% -------
% You may want to be able to select/exclude subsets of voxels and items without
% needed to make multiple copies of the data. By specifying filters, you can
% pre-specify these subsets and apply them programmatically.
% A filter is represented as a structure with 3 required fields: label,
% dimension, and filter. The label names the filter so that it can be easily
% referenced later, dimension encodes whether the filter applies to rows (1) or
% columns (2) of X. The filter is a binary vectory that represents the filter
% itself.
% Here, lets set up (totally arbitrarily) a ROI filter and a filter to exclude
% outliers.
z = [true(500,1);false(9500,1)];
FILTERS(1) = struct('label','ROI01', 'dimension', 2, 'filter', z);
z = [true(98,1);false(2,1)];
FILTERS(2) = struct('label','GoodRows', 'dimension', 1, 'filter', z);

% Coordinates
% -----------
% It is useful (and in the case of searchlight, essential) that the voxel
% coordinates be represented in the metadata structure. Like filters,
% coordinates are represented in a structure. The coordinate structure has 2
% required fields: 'orientation' and 'xyz'. The 'orientation' field is like the
% 'label' field in the filter structure and is used to look up particular
% coordinates. Labeling different coordinate spaces by 'orientation' is an AFNI
% convention. You don't have to use orientation codes like 'tlrc', 'orig', and
% 'mni', but that is my convention.
xyz = [(1:nVoxels)',ones(nVoxels,1),ones(nVoxels,1)];
COORDS(1) = struct('orientation','orig','xyz',xyz);
COORDS(2) = struct('orientation','tlrc','xyz',xyz);

% Put it all together
% -------------------
% The metadata object compiles these three items, along with a couple other
% bits, into a single structure. The metadata structure has several required
% fields: 'subject', 'targets', 'filters', 'coords', 'cvind', 'nrow', 'ncol',
% 'itemindex', and 'runindex'. The item index is used in case items are
% repeated and, for instance, may need to be averaged together. The run index
% is used to identify which block or scanner run each trial belongs to. There
% will be a metadata structure for each subject, compiled into a structured
% array. Although in the example below subjects 100 and 101 are the same aside
% from their subject numbers, in practice they could be given different
% information.
metadata(1).subject = 100;
metadata(1).targets = TARGETS;
metadata(1).filters = FILTERS;
metadata(1).coords = COORDS;
metadata(1).cvind = SCHEMES;
metadata(1).nrow = nItems;
metadata(1).ncol = nVoxels;
metadata(1).itemindex = 1:nItems;
metadata(1).runindex = ones(1,nItems);

metadata(2).subject = 101;
metadata(2).targets = TARGETS;
metadata(2).filters = FILTERS;
metadata(2).coords = COORDS;
metadata(2).cvind = SCHEMES;
metadata(2).nrow = nItems;
metadata(2).ncol = nVoxels;
metadata(2).itemindex = 1:nItems;
metadata(2).runindex = ones(1,nItems);

% Save the data to disk
% =====================
% Despite having data and metadata organized properly in memory, before working
% with WholeBrain_RSA we need to write the data to disk. The reason for this
% is that WholeBrain_RSA is not written to be used interactively, but rather
% to facilitate to use in headless, batch applications particularly on
% distributed computing systems. WholeBrain_RSA accepts paths to files on
% disk, as well as many other parameters.
% The data and metadata should be saved to a central location where it can be
% easily referenced.
% These files can be named whatever you like. You will be referencing them with
% explicit paths, and WholeBrain_RSA does not make any assumptions about them.
% The program does assume that the *variable* names are X and metadata, but
% this default can be overwritten with certain parameters to WholeBrain_RSA
% (data_var and metadata_var) if you prefer another convention.
subjects = [metadata.subject];
datadir = './shared';
if ~exist(datadir,'dir')
    mkdir(datadir);
end
Rotations = [30,60,90,120];
iRot = 0;
r = @(theta) [cos(theta), -sin(theta); sin(theta), cos(theta)];
for iSubj = 1:2
  s = subjects(iSubj);
  X = randn(nItems,nVoxels);
  % Rotate visual structure
  iRot = iRot + 1;
  theta = Rotations(iRot);
  Yv = Y.visual;
  n = size(Yv,2);
  R = eye(n);
  R(1:2, 1:2) = r(theta);
  Yv = Yv * R;
  % Rotate semantic structure
  iRot = iRot + 1;
  theta = Rotations(iRot);
  Ys = Y.semantic;
  n = size(Ys,2);
  R = eye(n);
  R(1:2, 1:2) = r(theta);
  Ys = Ys * R;
  % Add structure into X
  X(:,1:18) = X(:,1:18) + repmat(Yv,1,6);
  X(:,19:66) = X(:,19:66) + repmat(Ys,1,6);
  % Save data
  filename = sprintf('s%03d.mat', s);
  filepath = fullfile(datadir,filename);
  save(filepath, 'X');
end
save(fullfile(datadir,'metadata.mat'), 'metadata');

% Define a parameter file
% =======================
% WholeBrain_RSA, despite being written as a Matlab function, is a very ugly
% function. First of all, it does not return anything. All results are written
% to disk. Likewise, although it is possible to invoke WholeBrain_RSA from
% within a script or at the interactive terminal, it is designed to look for a
% parameter file if no arguments are provided. This is all makes
% WholeBrain_RSA a bit counter-intuitive. However, these design choices make
% much more sense when considered in a distributed computing environment.
% WholeBrain_RSA can be deployed to a system, along with a json file
% containing parameters, and it will parse the file and execute according to
% the instructions. It is designed to be executed with bare minimum interaction.
%
% Defining a parameter file is simple. See the documentation for a list of
% valid parameters. WholeBrain_RSA reads json (http://www.json.org/), which is
% a widely used text-based syntax for representing structured data.
%
%           **The file must be named params.json**
%
% To read and write json, you will need jsonlab
% (http://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files)
% which I have bundled with my code:
addpath('../dependencies/jsonlab/');

% Put the parameter file where you want to run the analysis. Paths can be
% relative with respect to where you execute WholeBrain_RSA, but in most cases
% it will probably make sense for them to be absolute. The following should
% translate into a valid json file for the purpose of this demo.
params = struct('regularization', 'GrOWL2', 'bias', false,'tau',0.2,...
    'lambda1', 0.4200556, 'lambda', 0.5863, 'LambdaSeq', 'inf',...
    'cvscheme', 1,'cvholdout', 1:10, 'finalholdout', 0,...
    'target', 'visual', 'sim_source', 'chamfer', 'sim_metric', 'chamfer',...
    'data', './shared/s100.mat', 'data_var', 'X',...
    'normalize', 'stdev', 'metadata', './shared/metadata.mat',...
    'metadata_varname', 'metadata', 'orientation', 'mni',...
    'filters', {{'ROI01','GoodRows'}}, 'SmallFootprint', false,...
    'debug', false, 'SaveResultsAs','json');

savejson('',params,'FileName','params.json','ForceRootName',false);

% Run WholeBrain_RSA
% ===================
% With data and metadata structured properly and saved to disk, and with a
% parameter file named params.json in a folder where you would like to execute
% the analysis and return results, all that remains is to boot up Matlab in the
% directory that contains 'params.json' and execute WholeBrain_RSA() at the
% command prompt. If you have compiled WholeBrain_RSA into an executable (as
% would be necessary on a distributed computing cluster), you can execute
% Wholebrain_MVPA directly from the command line. In either case, it will read
% the parameter file and begin analysis. When it completes you will find a
% results.mat (or results.json) file in the directory where WholeBrain_RSA was
% executed.
addpath('../src/')
WholeBrain_RSA()

% Compile Results
% ===============
% If you are using WholeBrain_RSA on a distributed computing cluster, you will
% quickly find that the volume of results is difficult to manage effectively. I
% have written some utility functions in Wholebrain_MVPA/util that attempt to
% facilitate common actions, like loading data from many jobs into a single
% matlab structure, writing tables of data, dumping coordinates of identified
% voxels, etc.
% Alternatively, you may find that your volume of data demands a database
% solution. Although the default is to return data in .mat files, which makes
% it easy to read back into matlab, results can also be output in json format
% which facilitates storing in a SQL or NoSQL database like MongoDB. Setting up
% such a database solution is far beyond the scope of this demo, but the squall
% project (github.com/ikinsella/squall) is a developing solution that utilizes
% MongoDB to great effect.
