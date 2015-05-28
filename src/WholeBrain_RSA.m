function WholeBrain_RSA()
  %@AVS = 'auditory', 'visual', 'semantic'
  %@lambda_in = scalar for group lasso and vector for group OWL (inversely proportional to sparsity)
  %@normalize = 1 for normalized feature vectors
  %@subject = '03' for subject 03

  % ----------------------Set parameters-----------------------------------------------
  jdat           = loadjson('params.json')
  DEBUG          = jdat.debug;
  Gtype          = jdat.Gtype;
  lambda_in      = jdat.lambda;
  normalize      = jdat.normalize
  BIAS           = jdat.bias
  AVS            = jdat.reptype;
  SimilarityType = jdat.SimilarityType;
  targets        = jdat.targets;
  datafile       = jdat.data;
  cvfile         = jdat.cvfile;
  cvholdout      = jdat.cvholdout;
  finalholdout   = jdat.finalholdout;
  metafile       = jdat.metadata;

  if isfield(jdat,'AdlasOpts')
    opts = jdat.AdlasOpts;
    % If values originated in a YAML file, and scientific notation is used, the
    % value may have been parsed as a string. Check and correct.
    if isfield(opts, 'tolInfeas')
        opts.tolInfeas = sscanf(opts.tolInfeas, '%e');
    end
    if isfield(opts, 'tolRelGap')
      if ~isnumeric(opts.tolRelGap)
        opts.tolRelGap = sscanf(opts.tolRelGap, '%e');
      end
    end
  else
    opts = struct();
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

  switch jdat.environment
  case 'condor'
    root = './';
    if isfield(jdat,'datadir')
      datadir = jdat.datadir
    else
      [ddir,dname,ext]  = fileparts(datafile);
      datadir           = fullfile(root,ddir);
      datafile          = strcat(dname,ext);
      [mdir,mname,ext]  = fileparts(metafile);
      datadir           = fullfile(root,mdir);
      metafile          = strcat(mname,ext);
      [cdir,cvname,ext] = fileparts(cvfile);
      cvdir             = fullfile(root,cdir);
      cvfile            = strcat(cvname,ext);
      assert(strcmp(mdir,ddir) && strcmp(cdir,mdir));
    end
    matfilename = 'results.mat';

  case 'chris'
    root = '/home/chris/MCW/WholeBrain_RSA';
    datadir = fullfile(root,'data');
    resultroot = fullfile(root,'results');
    if normalize == 1
      resultdir = fullfile(resultroot,'normalize',Gtype,AVS);
    else
      resultdir = fullfile(resultroot,'not_normalize',Gtype,AVS);
    end

    if DEBUG
      resultdir = fullfile(resultdir, 'DEBUG')
    end

    if ~exist(resultdir,'dir')
      mkdir(resultdir);
    end
    matfilename = fullfile(resultdir,sprintf('out_%d.mat',subjid));

  case 'urvashi'
    root = '/mnt/ws/home/uoswal/words/';
    datadir = fullfile(root,'data');
    resultroot = fullfile(root,'results');
    if normalize == 1
      resultdir = fullfile(resultroot,'normalize',Gtype,AVS);
    else
      resultdir = fullfile(resultroot,'not_normalize',Gtype,AVS);
    end

    if DEBUG
      resultdir = fullfile(resultdir, 'DEBUG')
    end

    if ~exist(resultdir,'dir')
      mkdir(resultdir);
    end
    matfilename = fullfile(resultdir,sprintf('out_%d.mat',subjid));

  otherwise
    error('Environment %s not implemented.', jdat.environment);

  end

  load(fullfile(datadir,metafile), 'metadata');
  [path,fname,ext] = fileparts(datafile);
  subjid = sscanf(fname, 's%d');
  subject = find(subjid == [metadata.subject]);


  warning off

  %Testing similarity based feature selection method
  % to fix/edit:
  % --extract similarity from distance matrix
  % --tune @lambda_try (and @R) parameters

  %% --------------------- Object-bank specific code------------------------------------

  switch targets
  case 'animals'
    ind = metadata(subject).TrueAnimals;

  case 'artifacts'
    ind = metadata(subject).TrueArtifacts;

  case 'all'
    ind = true(size(metadata(subject).TrueAnimals,1),1);

  otherwise
    error('Invalid target set: %s.', targets);

  end
  datapath = fullfile(datadir,datafile);
  fprintf('Loading data from  %s... ', datapath);
  load(datapath);

  % load input files and filter outliers more than 5 std dev away
  if isfield(jdat,'SanityCheckData')
    disp('PSYCH! This is a simulation.');
    switch jdat.SanityCheckData
    case 'shuffle'
      disp('Shuffling rows of MRI data!!!')
      X = X(randperm(size(X,1)),:);
    case 'random'
      disp('Generating totally random MRI data!!!')
      X = randn(size(X));
    case 'useshuffled'
      [path,fname,ext] = fileparts(datafile);
      fprintf('Using pre-shuffled MRI data!!! (for %s)\n', fname)
      rdatapath = fullfile(datadir, sprintf('%s_shuffle.mat',fname));
      load(rdatapath,'X')
    case 'userandom'
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
  ind = ind & reduxFilter.words';
  vox = allzero & reduxFilter.voxels;
  X = X(ind,vox);

  % Include voxel for bias
  if BIAS
    X = [X, ones(size(X,1),1)];
  end

  if isfield(jdat, 'cvfile')
    cvpath = fullfile(cvdir,cvfile);
    load(cvpath, 'CV');
    outlying_words = reduxFilter.words(ind);
    cvind = CV(outlying_words, jdat.cvscheme);
    finalholdout = cvind == jdat.finalholdout;
    X(finalholdout,:) = [];
    cvind(finalholdout) = [];
    cvind(cvind>jdat.finalholdout) = cvind(cvind>jdat.finalholdout) - 1;
    % Adjust the cv holdout index(es) down if they are higher than the final holdout.
    if ~isempty(cvholdout)
      cvholdout(cvholdout>jdat.finalholdout) = cvholdout(cvholdout>jdat.finalholdout) - 1;
    end

  else
    ncv   = jdat.ncv;
    [n,d] = size(X);
    cvind = mod(1:n, ncv)+1;
    cvind = cvind(randperm(n));
  end

  %% ----------------Visual, Audio or Semantic similarities and processing----------------

  switch AVS
  case 'auditory'
    error('AVS %s: Not implemented.', AVS)

  case 'visual'
    error('AVS %s: Not implemented.', AVS)

  case 'semantic'
    allSimStructs = load(fullfile(datadir,'semantic_model.mat'));
    sem = allSimStructs.(jdat.SimilarityType);
    S = sem(ind,ind); clear sem;
    S = S(~finalholdout, ~finalholdout);

  case 'orthographic'
    allSimStructs = load(fullfile(datadir,'orthography_model.mat'));
    orth = allSimStructs.(SimilarityType);
    S = orth(ind,ind); clear orth;
    S = S(~finalholdout, ~finalholdout);

  otherwise
    error('AVS %s not implemented. Should be aud, vis or sem.\n', AVS)

  end

  if isfield(jdat,'SanityCheckModel')
    switch jdat.SanityCheckModel
    case 'shuffle'
      disp('Shuffling Similarity Matrix!!!')
      shidx = randperm(size(S,1));
      S = S(shidx,shidx);
    case 'randcorr'
      disp('Generating totally random Similarity Matrix!!! (correlation)')
      x = randn(size(S,1),5);
      S = corr(x');
    case 'randinner'
      disp('Generating totally random Similarity Matrix!!! (inner product)')
      x = randn(size(S,1),5);
      S = x * x';
    case 'use_random_cor'
      disp('Using predefined random similarity matrix!!! (correlation)')
      [path,fname,ext] = fileparts(datafile);
      rdatapath = fullfile(datadir, 'random_model.mat');
      load(rdatapath,'cor');
      S = cor(ind,ind); clear cor;
      S = S(~finalholdout, ~finalholdout);
    case 'use_random_inner'
      disp('Using predefined random similarity matrix!!! (inner product)')
      [path,fname,ext] = fileparts(datafile);
      rdatapath = fullfile(datadir, 'random_model.mat');
      load(rdatapath,'inner');
      S = inner(ind,ind); clear inner;
      S = S(~finalholdout, ~finalholdout);
    case 'use_shuffled_cor'
      disp('Using pre-shuffled similarity matrix!!! (correlation)')
      [path,fname,ext] = fileparts(datafile);
      rdatapath = fullfile(datadir, 'shuffle_model.mat');
      load(rdatapath,'cor');
      S = cor(ind,ind); clear cor;
      S = S(~finalholdout, ~finalholdout);
    case 'use_shuffled_inner'
      disp('Using pre-shuffled similarity matrix!!! (inner product)')
      [path,fname,ext] = fileparts(datafile);
      rdatapath = fullfile(datadir, 'shuffle_model.mat');
      load(rdatapath,'inner');
      S = inner(ind,ind); clear inner;
      S = S(~finalholdout, ~finalholdout);
    case 'real'
      disp('Using the true similarity matrix, unaltered.')
    end
  end

  fprintf('Data loaded and processed.\n');

  %% ---------------------Setting algorithm parameters-------------------------
  [Uz, Sz, nz_rows, p1] = learn_similarity_encoding(S, X, lambda_in, cvind, cvholdout, normalize, Gtype, DEBUG, opts);

  fprintf('Saving stuff.....\n');

  summaryfile = fullfile(root,'summary.txt');

  save(matfilename,'p1','nz_rows','Sz','Uz');

  fprintf('Done!\n');
end
