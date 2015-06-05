function WholeBrain_RSA(varargin)
  p = inputParser;
  p.KeepUnmatched = true;
  % ----------------------Set parameters-----------------------------------------------
  addParameter(p , 'debug'            , false     , @islogicallike );
  addParameter(p , 'SmallFootprint'   , false     , @islogicallike );
  addParameter(p , 'Gtype'            , []        , @ischar        );
  addParameter(p , 'normalize'        , false     , @islogicallike );
  addParameter(p , 'bias'             , false     , @islogicallike );
  addParameter(p , 'simfile'          , []        , @ischar        );
  addParameter(p , 'simtype'          , []        , @ischar        );
  addParameter(p , 'filters'          , []                         );
  addParameter(p , 'data'             , []                         );
  addParameter(p , 'metadata'         , []        , @ischar        );
  addParameter(p , 'cvfile'           , []        , @ischar        );
  addParameter(p , 'cvscheme'         , []        , @isintegerlike );
  addParameter(p , 'cvholdout'        , []        , @isintegerlike );
  addParameter(p , 'finalholdout'     , []        , @isintegerlike );
  addParameter(p , 'tau'              , 0.2       , @isnumeric     );
  addParameter(p , 'lambda'           , []        , @isnumeric     );
  addParameter(p , 'lambda1'          , []        , @isnumeric     );
  addParameter(p , 'LambdaSeq'        , []        , @ischar        );
  addParameter(p , 'AdlasOpts'        , struct()  , @isstruct      );
  addParameter(p , 'environment'      , 'condor'  , @ischar        );
  addParameter(p , 'SanityCheckData'  , []        , @ischar        );
  addParameter(p , 'SanityCheckModel' , []        , @ischar        );

  if nargin > 0
    parse(p, varargin{:});
  else
    jdat = loadjson('params.json');
    fields = fieldnames(jdat);
    jcell = [fields'; struct2cell(jdat)']
    parse(p, jcell{:});
  end

  % private function.
  assertRequiredParameters(p.Results);

  DEBUG            = p.Results.debug;
  SmallFootprint   = p.Results.SmallFootprint;
  Gtype            = p.Results.Gtype;
  normalize        = p.Results.normalize;
  BIAS             = p.Results.bias;
  simfile          = p.Results.simfile;
  simtype          = p.Results.simtype;
  filters          = p.Results.filters;
  datafile         = p.Results.data;
  cvfile           = p.Results.cvfile;
  cvscheme         = p.Results.cvscheme;
  cvholdout        = p.Results.cvholdout;
  finalholdoutInd  = p.Results.finalholdout;
  metafile         = p.Results.metadata;
  tau              = p.Results.tau;
  lambda           = p.Results.lambda;
  lambda1          = p.Results.lambda1;
  LambdaSeq        = p.Results.LambdaSeq;
  opts             = p.Results.AdlasOpts;
  environment      = p.Results.environment;
  SanityCheckData  = p.Results.SanityCheckData;
  SanityCheckModel = p.Results.SanityCheckModel;

  % Check that the correct parameters are passed, given the desired algorithm
  [lam, lam1, lamSeq] = verifyLambdaSetup(Gtype, lambda, lambda1, LambdaSeq);

  % If values originated in a YAML file, and scientific notation is used, the
  % value may have been parsed as a string. Check and correct.
  if isfield(opts, 'tolInfeas')
    if ischar(opts.tolInfeas)
      opts.tolInfeas = sscanf(opts.tolInfeas, '%e');
    end
  end
  if isfield(opts, 'tolRelGap')
    if ischar(opts.tolRelGap)
      opts.tolRelGap = sscanf(opts.tolRelGap, '%e');
    end
  end

  if iscell(datafile)
    if length(datafile) == 1
      datafile = datafile{1};
    end
  end
  if iscell(metafile)
    if length(metafile) == 1
      metafile = metafile{1};
    end
  end

  switch environment
  case 'condor'
    root = './';
    datadir = root;
    matfilename = 'results.mat';
    infofilename = 'info.mat';
  case 'chris'
    root = './';
    datadir = fullfile(root,'data');
    matfilename = 'results.mat';
    infofilename = 'info.mat';

  otherwise
    error('Environment %s not implemented.', environment);

  end

  load(fullfile(datadir,metafile), 'metadata');
  [path,fname,ext] = fileparts(datafile);
  subjid = sscanf(fname, 's%d');
  if ~isempty(subjid)
    subject = find(subjid == [metadata.subject]);
  end

  warning off

  %Testing similarity based feature selection method
  % to fix/edit:
  % --extract similarity from distance matrix
  % --tune @lambda_try (and @R) parameters

  %% --------------------- Object-bank specific code------------------------------------
  datapath = fullfile(datadir,datafile);
  fprintf('Loading data from  %s...\n', datapath);
  load(datapath, 'X');

  if ~isempty(filters)
    if iscell(filters)
        n = length(filters);
        key = filters{1};
    else
        n = 1;
        key = filters;
    end
    z = metadata.(key);
    filter = z;
    if n > 1
      for i = 2:n
        key = filters{i};
        z = metadata.(key);
        filter(z) = true;
      end
    end
  else
    filter = true(size(X,1),1);
  end

  % load input files and filter outliers more than 5 std dev away
  if ~isempty(SanityCheckData)
    disp('PSYCH! This is a simulation.');
    switch SanityCheckData
    case 'shuffle'
      disp('Shuffling rows of MRI data!!!')
      X = X(randperm(size(X,1)),:);
    case 'random'
      disp('Generating totally random MRI data!!!')
      X = randn(size(X));
    case 'use_shuffled'
      [path,fname,ext] = fileparts(datafile);
      fprintf('Using pre-shuffled MRI data!!! (for %s)\n', fname)
      rdatapath = fullfile(datadir, sprintf('%s_shuffle.mat',fname));
      load(rdatapath,'X')
    case 'use_random'
      [path,fname,ext] = fileparts(datafile);
      fprintf('Using predefined random MRI data!!! (for %s)\n', fname)
      rdatapath = fullfile(datadir, sprintf('%s_random.mat',fname));
      load(rdatapath,'X')
    case 'real'
      disp('Using the true data, unaltered.');
    end
  end
  allzero = any(X); % Identify columns with data
  [~, reduxFilter] = removeOutliers(X);
  % Note: reduxFilter field names are reversed...
  filter = filter & reduxFilter.words';
  vox = allzero & reduxFilter.voxels;
  X = X(filter,vox);

  % Include voxel for bias
  fprintf('%-28s', 'Including Bias Unit:');
  msg = 'NO';
  if BIAS
    msg = 'YES';
    X = [X, ones(size(X,1),1)];
  end
  fprintf('[%3s]\n', msg);

  % Normalize columns of X
  fprintf('%-28s', 'Normalizing columns of X:');
  msg = 'NO';
  if normalize
    msg = 'YES';
  end
  fprintf('[%3s]\n', msg);

  cvpath = fullfile(datadir,cvfile);
  load(cvpath, 'CV');
  outlying_words = reduxFilter.words(filter);
  cvind = CV(outlying_words, cvscheme);
  finalholdout = cvind == finalholdoutInd;
  X(finalholdout,:) = [];
  cvind(finalholdout) = [];
  if finalholdoutInd > 0
    cvind(cvind>finalholdoutInd) = cvind(cvind>finalholdoutInd) - 1;
    % Adjust the cv holdout index(es) down if they are higher than the final holdout.
    if ~isempty(cvholdout)
      cvholdout(cvholdout>finalholdoutInd) = cvholdout(cvholdout>finalholdoutInd) - 1;
    end
  end

  %% ----------------Visual, Audio or Semantic similarities and processing----------------

  if ~isempty(SanityCheckModel)
    switch SanityCheckModel
    case 'shuffle'
      disp('Shuffling Similarity Matrix!!!')
      simpath = fullfile(datadir,simfile);
      allSimStructs = load(simpath);
      S = allSimStructs.(simtype);
      shidx = randperm(size(S,1));
      S = S(shidx,shidx);

    case 'random'
      fprintf('Generating totally random Similarity Matrix!!! ')
      simpath = fullfile(datadir,simfile);
      allSimStructs = load(simpath);
      S = allSimStructs.(simtype);
      x = randn(size(S,1),5);

      switch simtype
      case 'inner'
        fprintf('(inner product)\n');
        S = x * x';
      case 'cor'
        fprintf('(correlation)\n');
        S = corr(x');
      case 'cosine'
        fprintf('(cosine)\n');
        S = cosinesimilarity(x);
      case 'earthmover'
        fprintf('(cosine---EMD not implemented)\n');
        S = cosinesimilarity(x);
      end

    case 'use_random'
      fprintf('Using predefined random similarity matrix!!! ')
      switch simtype
      case 'inner'
        fprintf('(inner product)\n');
      case 'cor'
        fprintf('(correlation)\n');
      case 'cosine'
        fprintf('(cosine)\n');
      case 'earthmover'
        fprintf('(cosine---EMD not implemented)\n');
      end
      [~,fname,~] = fileparts(simfile);
      simfile = sprintf('%s_random.mat',fname);
      simpath = fullfile(datadir, simfile);
      allSimStructs = load(simpath);
      S = allSimStructs.(simtype);

    case 'use_shuffled'
      fprintf('Using pre-shuffled similarity matrix!!! ')
      switch simtype
      case 'inner'
        fprintf('(inner product)\n');
      case 'cor'
        fprintf('(correlation)\n');
      case 'cosine'
        fprintf('(cosine)\n');
      case 'earthmover'
        fprintf('(Earth Movers Distance)\n');
      end
      [~,fname,~] = fileparts(simfile);
      simfile = sprintf('%s_shuffled.mat',fname);
      simpath = fullfile(datadir, simfile);
      allSimStructs = load(simpath);
      S = allSimStructs.(simtype);

    case 'real'
      disp('Using the true similarity matrix, unaltered.')
      simpath = fullfile(datadir,simfile);
      allSimStructs = load(simpath);
      S = allSimStructs.(simtype);
    end
  else
    simpath = fullfile(datadir,simfile);
    allSimStructs = load(simpath);
    S = allSimStructs.(simtype);
  end
  S = S(filter,filter); clear allSimStructs;
  S = S(~finalholdout, ~finalholdout);

  fprintf('Data loaded and processed.\n');

  %% ---------------------Setting algorithm parameters-------------------------
  switch Gtype
  case 'L1L2'
    [results,info] = learn_similarity_encoding(S, X, Gtype, ...
                      'tau'            , tau            , ...
                      'lambda1'        , lambda1        , ...
                      'cvind'          , cvind          , ...
                      'cvholdout'      , cvholdout      , ...
                      'normalize'      , normalize      , ...
                      'DEBUG'          , DEBUG          , ...
                      'SmallFootprint' , SmallFootprint , ...
                      'AdlasOpts'      , opts); %#ok<ASGLU>

  case 'grOWL'
    [results,info] = learn_similarity_encoding(S, X, Gtype, ...
                      'tau'            , tau            , ...
                      'lambda'         , lambda         , ...
                      'LambdaSeq'      , LambdaSeq      , ...
                      'cvind'          , cvind          , ...
                      'cvholdout'      , cvholdout      , ...
                      'normalize'      , normalize      , ...
                      'DEBUG'          , DEBUG          , ...
                      'SmallFootprint' , SmallFootprint , ...
                      'AdlasOpts'      , opts); %#ok<ASGLU>

  case 'grOWL2'
    [results,info] = learn_similarity_encoding(S, X, Gtype, ...
                      'tau'            , tau            , ...
                      'lambda'         , lambda         , ...
                      'lambda1'        , lambda1        , ...
                      'LambdaSeq'      , LambdaSeq      , ...
                      'cvind'          , cvind          , ...
                      'cvholdout'      , cvholdout      , ...
                      'normalize'      , normalize      , ...
                      'DEBUG'          , DEBUG          , ...
                      'SmallFootprint' , SmallFootprint , ...
                      'AdlasOpts'      , opts); %#ok<ASGLU>
  end

  fprintf('Saving stuff.....\n');

  save(matfilename,'-struct','results');
  save(infofilename,'-struct','info');

  fprintf('Done!\n');
end

function [lam, lam1, lamSeq] = verifyLambdaSetup(Gtype, lambda, lambda1, LambdaSeq);
% Each algorithm requires different lambda configurations. This private
% function ensures that everything has been properly specified.
  switch Gtype
  case 'L1L2'
    if ~isempty(lambda)
      warning('Group Lasso does not use the lambda parameter. It is being ignored.');
    end
    assert(~isempty(lambda1)   , 'Group Lasso requires lambda1.');
    lam    = [];
    lam1   = lambda1;
    lamSeq = [];

  case 'grOWL'
    if ~isempty(lambda1)
      warning('grOWL does not use the lambda1 parameter. It is being ignored.');
    end
    assert(~isempty(lambda)    , 'Group Lasso requires lambda.');
    assert(~isempty(LambdaSeq) , 'A LambdaSeq type (linear or exponential) must be set when using grOWL*.');
    lam    = lambda;
    lam1   = [];
    lamSeq = LambdaSeq;

  case 'grOWL2'
    assert(~isempty(lambda)    , 'grOWL2 Lasso requires lambda.');
    assert(~isempty(lambda1)   , 'grOWL2 requires lambda1.');
    assert(~isempty(LambdaSeq) , 'A LambdaSeq type (linear or exponential) must be set when using grOWL*.');
    lam    = lambda;
    lam1   = lambda1;
    lamSeq = LambdaSeq;
  end
end

function assertRequiredParameters(params)
  required = {'Gtype'    , 'simfile' , 'simtype'   , 'data'     , ...
              'metadata' , 'cvfile'  , 'cvscheme'  ,'cvholdout' , 'finalholdout'};
  N = length(required);
  for i = 1:N
    req = required{i};
    assert(isfield(params,req), '%s must exist in params structure! Exiting.');
    assert(~isempty(params.(req)), '%s must be set. Exiting.');
  end
end

function b = islogicallike(x)
  b = any(x == [1,0]);
end

function b = isintegerlike(x)
  b = mod(x,1) == 0;
end
