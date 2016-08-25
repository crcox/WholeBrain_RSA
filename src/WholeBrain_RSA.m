function WholeBrain_RSA(varargin)
  p = inputParser;
  p.KeepUnmatched = false;
  % ----------------------Set parameters-----------------------------------------------
  addParameter(p , 'debug'            , false     , @islogicallike );
  addParameter(p , 'RandomSeed'       , 0                          );
  addParameter(p , 'PermutationTest'  , false     , @islogicallike );
  addParameter(p , 'SmallFootprint'   , false     , @islogicallike );
  addParameter(p , 'regularization'   , []        , @ischar        );
  addParameter(p , 'normalize'        , false                      );
  addParameter(p , 'bias'             , false     , @islogicallike );
  addParameter(p , 'target'           , []        , @ischar        );
  addParameter(p , 'target_type'      , []        , @ischar        );
  addParameter(p , 'sim_source'       , []        , @ischar        );
  addParameter(p , 'sim_metric'       , []        , @ischar        );
  addParameter(p , 'filters'          , []                         );
  addParameter(p , 'data'             , []                         );
  addParameter(p , 'data_varname'     , []                         );
  addParameter(p , 'metadata'         , []        , @ischar        );
  addParameter(p , 'metadata_varname' , []        , @ischar        );
  addParameter(p , 'finalholdout'     , 0         , @isintegerlike );
  addParameter(p , 'cvscheme'         , []        , @isnumeric     );
  addParameter(p , 'cvholdout'        , []        , @isnumeric     );
  addParameter(p , 'orientation'      , []        , @ischar        );
  addParameter(p , 'tau'              , 0.2       , @isnumeric     );
  addParameter(p , 'lambda'           , []        , @isnumeric     );
  addParameter(p , 'lambda1'          , []        , @isnumeric     );
  addParameter(p , 'LambdaSeq'        , []        , @ischar        );
  addParameter(p , 'AdlasOpts'        , struct()  , @isstruct      );
  addParameter(p , 'SanityCheckData'  , []        , @ischar        );
  addParameter(p , 'SanityCheckModel' , []        , @ischar        );
  addParameter(p , 'SaveResultsAs'  , 'mat'       , @isMatOrJSON   );
  addParameter(p , 'subject_id_fmt' , '%d'        , @ischar        );
  % --- searchlight specific --- %
  addParameter(p , 'searchlight'      , 0        , @islogicallike );
  addParameter(p , 'slShape'          , ''        , @ischar        );
  addParameter(p , 'slSim_Measure'    , ''        , @ischar        );
  addParameter(p , 'slRadius'         , []        , @isnumeric     );
  addParameter(p , 'slPermutationType', ''        , @ischar        );
  addParameter(p , 'slPermutations'   , 0         , @isscalar      );
  % Parameters below this line are unused in the analysis, may exist in the
  % parameter file because other progams use them.
  addParameter(p , 'COPY'             , []                         );
  addParameter(p , 'URLS'             , []                         );
  addParameter(p , 'executable'       , []                         );
  addParameter(p , 'wrapper'          , []                         );

  if nargin > 0
    parse(p, varargin{:});
  else
    jdat = loadjson('params.json');
    fields = fieldnames(jdat);
    jcell = [fields'; struct2cell(jdat)'];
    parse(p, jcell{:});
  end


  % private function.
  assertRequiredParameters(p.Results);

  DEBUG            = p.Results.debug;
  PermutationTest  = p.Results.PermutationTest;
  SmallFootprint   = p.Results.SmallFootprint;
  RandomSeed       = p.Results.RandomSeed;
  regularization   = p.Results.regularization;
  normalize        = p.Results.normalize;
  BIAS             = p.Results.bias;
  target_label     = p.Results.target;
  target_type      = p.Results.target_type;
  sim_source       = p.Results.sim_source;
  sim_metric       = p.Results.sim_metric;
  filter_labels    = p.Results.filters;
  datafile         = p.Results.data;
  data_varname     = p.Results.data_varname;
  cvscheme         = p.Results.cvscheme;
  cvholdout        = p.Results.cvholdout;
  finalholdoutInd  = p.Results.finalholdout;
  orientation      = p.Results.orientation;
  metafile         = p.Results.metadata;
  metadata_varname = p.Results.metadata_varname;
  tau              = p.Results.tau;
  lambda           = p.Results.lambda;
  lambda1          = p.Results.lambda1;
  LambdaSeq        = p.Results.LambdaSeq;
  opts             = p.Results.AdlasOpts;
  SanityCheckData  = p.Results.SanityCheckData;
  SanityCheckModel = p.Results.SanityCheckModel;
  SaveResultsAs    = p.Results.SaveResultsAs;
  FMT_subjid       = p.Results.subject_id_fmt;
  % --- searchlight specific --- %
  SEARCHLIGHT   = p.Results.searchlight;
  slSim_Measure = p.Results.slSim_Measure;
  slPermutationType = p.Results.slPermutationType;
  slPermutationCount = p.Results.slPermutations;
  slShape = p.Results.slShape;
  slRadius = p.Results.slRadius;

  rng(RandomSeed);

  % Check that the correct parameters are passed, given the desired regularization
  if ~isempty(regularization)
    [lambda, lambda1, LambdaSeq] = verifyLambdaSetup(regularization, lambda, lambda1, LambdaSeq);
  end
  if SEARCHLIGHT && ~strcmpi(slSim_Measure,'nrsa')
    assert(~isempty(slPermutationType));
    assert(~isempty(slPermutationCount));
  end

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

  % If cell array with one element, unpack element from cell.
  datafile = uncell(datafile);
  metafile = uncell(metafile);

  %% Load metadata
  StagingContainer = load(metafile, metadata_varname);
  metadata = StagingContainer.(metadata_varname); clear StagingContainer;
  N = length(metadata);
  n = [metadata.nrow];
  d = [metadata.ncol];

  %% Compile filters
  rowfilter  = cell(N,1);
  colfilter  = cell(N,1);
  for i = 1:N
    if isempty(filter_labels)
      rowfilter{i} = true(1,n(i));
      colfilter{i} = true(1,d(i));
    else
      [rowfilter{i},colfilter{i}] = composeFilters(metadata(i).filters, filter_labels);
    end
  end

  if ~isempty(SEARCHLIGHT) && SEARCHLIGHT && strcmpi(slSim_Measure,'nrsa') && finalholdoutInd > 0
    %% Select targets
    Sall = selectTargets(metadata, target_type, target_label, sim_source, sim_metric, rowfilter);

    %% Load data
    [Xall,subjix] = loadData(datafile, data_varname, rowfilter, colfilter, metadata);
    if iscell(Xall) && numel(Xall) == 1
      Xall = Xall{1};
    end
    Sall = Sall{subjix};
  end

  %% Load CV indexes, and identify the final holdout set.
  % N.B. the final holdout set is excluded from the rowfilter.
  cvind = cell(1,N);
  cvindAll = cell(1,N);
  for i = 1:N
    % Add the final holdout set to the rowfilter, so we don't even load
    % those data.
    cvindAll{i} = metadata(i).cvind(:,cvscheme);
    finalholdout = cvindAll{i} == finalholdoutInd;
    % Remove the final holdout set from the cvind, to match.
    rowfilter{i} = forceRowVec(rowfilter{i}) & forceRowVec(~finalholdout);
    cvind{i} = cvindAll{i}(rowfilter{i});
  end

  %% Select targets
  fprintf('\n');
  fprintf('Loading similarity structure\n');
  fprintf('----------------------------\n');
  fprintf('%12s: %s\n', 'target_label', target_label);
  fprintf('%12s: %s\n', 'type', target_type);
  fprintf('%12s: %s\n', 'sim_source', sim_source);
  fprintf('%12s: %s\n', 'sim_metric', sim_metric);
  fprintf('\n');
  S = selectTargets(metadata, target_type, target_label, sim_source, sim_metric, rowfilter);

  %% Load data
  [X,subjix] = loadData(datafile, data_varname, rowfilter, colfilter, metadata, FMT_subjid);
  if iscell(X) && numel(X) == 1
    X = X{1};
  end
  S = S{subjix};
  metadata   = metadata(subjix);
  rowfilter  = rowfilter{subjix};
  colfilter  = colfilter{subjix};
  cvind      = cvind{subjix};
  cvindAll   = cvindAll{subjix};
  z = strcmp({metadata.coords.orientation}, orientation);
  COORDS = metadata.coords(z);
  COORDS_FIELDS = fieldnames(COORDS);
  for i = 1:numel(COORDS_FIELDS)
    cfield = COORDS_FIELDS{i};
    switch cfield
    case 'ind'
      COORDS.ind = COORDS.ind(colfilter);
    case 'ijk'
      COORDS.ijk = COORDS.ijk(colfilter,:);
    case 'xyz'
      COORDS.xyz = COORDS.xyz(colfilter,:);
    end
  end
  xyz = COORDS.xyz;
  fprintf('Initial dimensions: (%d,%d)\n', size(X,1), size(X,2));
  fprintf('Filtered dimensions: (%d,%d)\n', size(X,1), size(X,2));

  %% Include voxel for bias
  fprintf('%-28s', 'Including Bias Unit:');
  msg = 'NO';
  if BIAS
    msg = 'YES';
    X = [X, ones(size(X,1),1)];
  end
  fprintf('[%3s]\n', msg);

  %% Normalize columns of X
  % NB The normalization happens later.
  fprintf('%-28s', 'Normalizing columns of X:');
  msg = 'NO';
  if normalize
    msg = 'YES';
  end
  fprintf('[%3s]\n', msg);

  fprintf('Data loaded and processed.\n');

  %% ---------------------Setting regularization parameters-------------------------
  if SEARCHLIGHT
    X = uncell(X);
    S = uncell(S);
    cvind = uncell(cvind);
    cvset = unique(cvind);
    colfilter = uncell(colfilter);

    % create a 3D binary mask
    [mask,dxyz] = coordsTo3dMask(xyz);

    % Translate slradius (in mm) to sl voxels
    % N.B. Because voxels need not be symmetric cubes, but Seachmight will
    % generate symmetric spheres from a single radius parameter, we need to
    % select one value of the three that will be produced in this step. I am
    % arbitrarily choosing the max, to err on the side of being inclusive.
    slradius_ijk = max(round(slRadius ./ dxyz));

    % create the "meta" neighbourhood structure
    meta = createMetaFromMask(mask, 'radius', slradius_ijk);
    labels = metadata.itemindex(rowfilter);
    labelsRun = metadata.runindex(rowfilter);

    results.similarity_measure = slSim_Measure;
    if strcmpi('nrsa',slSim_Measure)
      % Define results structure
      results.Uz = [];
      results.Cz = [];
      results.Sz = [];
      results.nz_rows =  [];
      results.target_label = target_label;
      results.subject =  [];
      results.cvholdout = [];
      results.finalholdout = [];
      results.lambda = [];
      results.lambda1 = [];
      results.LambdaSeq = [];
      results.regularization = [];
      results.bias = [];
      results.normalize = [];
      results.nzv = [];
%      results.p1      =  [];
%      results.p2      =  [];
%      results.cor1    =  [];
%      results.cor2    =  [];
%      results.p1t     =  [];
%      results.p2t     =  [];
%      results.cor1t   =  [];
%      results.cor2t   =  [];
      results.coords  = [];
      results.structureScoreMap = zeros(1, size(meta.voxelsToNeighbours,1));
      results.structurePvalueMap = zeros(1, size(meta.voxelsToNeighbours,1));
      results.err1    =  zeros(1, size(meta.voxelsToNeighbours,1));
      results.err2    =  zeros(1, size(meta.voxelsToNeighbours,1));
      results.iter    =  [];

      % Preallocate
      if isempty(lambda); nlam = 1; else nlam = numel(lamba); end
      if isempty(lambda); nlam1 = 1; else nlam1 = numel(lambda1); end
      results(numel(cvset)*nlam*nlam1).Uz = [];

      for iVolume = 1:size(meta.voxelsToNeighbours,1)
        sl = meta.voxelsToNeighbours(iVolume,1:meta.numberOfNeighbours(iVolume));
        switch upper(regularization)
        case 'L1L2'
          [lambda1, err_L1L2] = fminbnd(@(x) optimizeGroupLasso(S,X(:,sl),tau,cvind,cvholdout,normalize,PermutationTest,x), 0, 32);
          if finalholdout > 0
            [tmpr,info] = learn_similarity_encoding(Sall, Xall(:,sl), regularization, target_type,...
                              'tau'            , tau            , ...
                              'lambda1'        , lambda1        , ...
                              'cvind'          , cvindAll       , ...
                              'cvholdout'      , finalholdoutInd, ...
                              'normalize'      , normalize      , ...
                              'DEBUG'          , DEBUG          , ...
                              'PermutationTest', PermutationTest, ...
                              'SmallFootprint' , SmallFootprint , ...
                              'AdlasOpts'      , opts); %#ok<ASGLU>
          else
            [tmpr,info] = learn_similarity_encoding(S, X(:,sl), regularization, target_type,...
                              'tau'            , tau            , ...
                              'lambda1'        , lambda1        , ...
                              'cvind'          , cvind          , ...
                              'cvholdout'      , cvholdout      , ...
                              'normalize'      , normalize      , ...
                              'DEBUG'          , DEBUG          , ...
                              'PermutationTest', PermutationTest, ...
                              'SmallFootprint' , SmallFootprint , ...
                              'AdlasOpts'      , opts); %#ok<ASGLU>
          end

        case 'GROWL'
          [tmpr,info] = learn_similarity_encoding(S, X(:,sl), regularization, target_type,...
                            'tau'            , tau            , ...
                            'lambda'         , lambda         , ...
                            'LambdaSeq'      , LambdaSeq      , ...
                            'cvind'          , cvind          , ...
                            'cvholdout'      , cvholdout      , ...
                            'normalize'      , normalize      , ...
                            'DEBUG'          , DEBUG          , ...
                            'PermutationTest', PermutationTest, ...
                            'SmallFootprint' , SmallFootprint , ...
                            'AdlasOpts'      , opts); %#ok<ASGLU>

        case 'GROWL2'
          [tmpr,info] = learn_similarity_encoding(S, X(:,sl), regularization, target_type,...
                            'tau'            , tau            , ...
                            'lambda'         , lambda         , ...
                            'lambda1'        , lambda1        , ...
                            'LambdaSeq'      , LambdaSeq      , ...
                            'cvind'          , cvind          , ...
                            'cvholdout'      , cvholdout      , ...
                            'normalize'      , normalize      , ...
                            'DEBUG'          , DEBUG          , ...
                            'PermutationTest', PermutationTest, ...
                            'SmallFootprint' , SmallFootprint , ...
                            'AdlasOpts'      , opts); %#ok<ASGLU>
        end
        for iResult = 1:numel(tmpr)
          results(iResult).err1(iVolume) = tmpr(iResult).err1;
          results(iResult).err2(iVolume) = tmpr(iResult).err2;
          results(iResult).structureScoreMap(iVolume) = tmpr(iResult).structureScoreMap;
        end
      end
    else
      fprintf('PermutationTest: %d\n', PermutationTest);
      if PermutationTest
        for ic = unique(cvind)'
          fprintf('Permuting CV %d...\n', ic);
          s = S(cvind==ic, cvind==ic);
          n = size(s,1);
          permix = randperm(n);
          S(cvind==ic, cvind==ic) = S(permix, permix);
        end
      end

      [structureScoreMap] = computeSimilarityStructureMap(...
        slSim_Measure,...
        X,labels,...
        X,labels,...
        'meta',meta,'similarityStructure',S,...
        'permutationTest',slPermutationType, slPermutationCount,...
        'groupLabels',labelsRun,labelsRun);

      results.structureScoreMap = structureScoreMap;
      results.RandomSeed = RandomSeed;
    end

    for iResult = 1:numel(results)
      results(iResult).coords = COORDS;
    end
  else
    switch upper(regularization)
    case 'L1L2'
      [results,info] = learn_similarity_encoding(S, X, regularization, target_type,...
                        'tau'            , tau            , ...
                        'lambda1'        , lambda1        , ...
                        'cvind'          , cvind          , ...
                        'cvholdout'      , cvholdout      , ...
                        'normalize'      , normalize      , ...
                        'DEBUG'          , DEBUG          , ...
                        'PermutationTest', PermutationTest, ...
                        'SmallFootprint' , SmallFootprint , ...
                        'AdlasOpts'      , opts); %#ok<ASGLU>

    case 'GROWL'
      [results,info] = learn_similarity_encoding(S, X, regularization, target_type,...
                        'tau'            , tau            , ...
                        'lambda'         , lambda         , ...
                        'lambda1'        , lambda1        , ...
                        'LambdaSeq'      , LambdaSeq      , ...
                        'cvind'          , cvind          , ...
                        'cvholdout'      , cvholdout      , ...
                        'normalize'      , normalize      , ...
                        'DEBUG'          , DEBUG          , ...
                        'PermutationTest', PermutationTest, ...
                        'SmallFootprint' , SmallFootprint , ...
                        'AdlasOpts'      , opts); %#ok<ASGLU>

    case 'GROWL2'
      [results,info] = learn_similarity_encoding(S, X, regularization, target_type,...
                        'tau'            , tau            , ...
                        'lambda'         , lambda         , ...
                        'lambda1'        , lambda1        , ...
                        'LambdaSeq'      , LambdaSeq      , ...
                        'cvind'          , cvind          , ...
                        'cvholdout'      , cvholdout      , ...
                        'normalize'      , normalize      , ...
                        'DEBUG'          , DEBUG          , ...
                        'PermutationTest', PermutationTest, ...
                        'SmallFootprint' , SmallFootprint , ...
                        'AdlasOpts'      , opts); %#ok<ASGLU>
    end
    for iResult = 1:numel(results)
      results(iResult).coords = COORDS;
      for i = 1:numel(COORDS_FIELDS)
        cfield = COORDS_FIELDS{i};
        switch cfield
        case 'ind'
          tmpind = COORDS.ind(results(iResult).Uix);
          results(iResult).coords.ind = tmpind(:)';
        case 'ijk'
          results(iResult).coords.ijk = COORDS.ijk(results(iResult).Uix,:);
        case 'xyz'
          results(iResult).coords.xyz = COORDS.xyz(results(iResult).Uix,:);
        end
      end
    end
  end

  fprintf('Saving stuff.....\n');

  [results.subject] = deal(subjix);
  [results.finalholdout] = deal(finalholdoutInd);

  %% Save results
  rinfo = whos('results');
  switch SaveResultsAs
    case 'mat'
      if rinfo.bytes > 2e+9 % 2 GB
        save('results.mat','results','-v7.3');
      else
        save('results.mat','results');
      end
    case 'json'
      if rinfo.bytes > 16e+6 % 16 MB
        disp('WARNING: Results structure too large to save as JSON (excedes MongoDB 16MB limit). Saving as .mat...')
        if rinfo.bytes > 2e+9 % 2 GB
          save('results.mat','results','-v7.3');
        else
          save('results.mat','results');
        end
      else
        savejson('',results,'FileName','results.json','ForceRootName',false);
      end
  end

  fprintf('Done!\n');
end

function [lam, lam1, lamSeq] = verifyLambdaSetup(regularization, lambda, lambda1, LambdaSeq)
% Each regularization requires different lambda configurations. This private
% function ensures that everything has been properly specified.
  switch upper(regularization)
  case 'NONE'
    if ~isempty(lambda) || ~isempty(lambda1)
      warning('Regularization was set to none, but lambda values were provided. They will be ignored.')
    end
    lam    = [];
    lam1   = [];
    lamSeq = [];

  case 'L1L2'
    if ~isempty(lambda)
      warning('Group Lasso does not use the lambda parameter. It is being ignored.');
    end
    assert(~isempty(lambda1)   , 'Group Lasso requires lambda1.');
    lam    = [];
    lam1   = lambda1;
    lamSeq = [];

  case 'GROWL'
    assert(~isempty(lambda) && ~isnan(lambda), 'grOWL requires lambda.');
    assert(~isempty(lambda1) && ~isnan(lambda1), 'grOWL requires lambda1.');
    assert(~isempty(LambdaSeq), 'A LambdaSeq type (linear or exponential) must be set when using grOWL*.');
    lam    = lambda;
    lam1   = lambda1;
    lamSeq = LambdaSeq;

  case 'GROWL2'
    assert(~isempty(lambda)    , 'grOWL2 requires lambda.');
    assert(~isempty(lambda1)   , 'grOWL2 requires lambda1.');
    assert(~isempty(LambdaSeq) , 'A LambdaSeq type (linear or exponential) must be set when using grOWL*.');
    lam    = lambda;
    lam1   = lambda1;
    lamSeq = LambdaSeq;
  end
end

function assertRequiredParameters(params)
  required = {'target','sim_metric','sim_source','data', ...
              'metadata','cvscheme','cvholdout','finalholdout','orientation'};
  N = length(required);
  for i = 1:N
    req = required{i};
    assert(isfield(params,req), '%s must exist in params structure! Exiting.',req);
    assert(~isempty(params.(req)), '%s must be set. Exiting.',req);
  end
end

function b = islogicallike(x)
  b = any(x == [1,0]);
end

function b = isintegerlike(x)
  b = mod(x,1) == 0;
end

function b = isMatOrJSON(x)
    b = any(strcmpi(x, {'mat','json'}));
end
