function WholeBrain_RSA()

  % ----------------------Set parameters-----------------------------------------------
  jdat         = loadjson('params.json')
  DEBUG        = jdat.debug;
  Gtype        = jdat.Gtype;
  lambda_in    = jdat.lambda;
  normalize    = jdat.normalize;
  BIAS         = jdat.bias;
  simfile      = jdat.simfile;
  simtype      = jdat.simtype;
  filters      = jdat.filters;
  datafile     = jdat.data;
  cvfile       = jdat.cvfile;
  cvholdout    = jdat.cvholdout;
  finalholdout = jdat.finalholdout;
  if isfield(jdat,'metadata')
    metafile       = jdat.metadata;
  end

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
  if exist('metafile','var')
    if iscell(metafile)
      if length(metafile) == 1
        metafile = metafile{1};
      end
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
      [cdir,cvname,ext] = fileparts(cvfile);
      cvdir             = fullfile(root,cdir);
      cvfile            = strcat(cvname,ext);
      [mdir,mname,ext]  = fileparts(metafile);
      datadir           = fullfile(root,mdir);
      metafile          = strcat(mname,ext);
      assert(strcmp(mdir,ddir) && strcmp(cdir,mdir));
    end
    matfilename = 'results.mat';

  otherwise
    error('Environment %s not implemented.', jdat.environment);

  end

  if exist('metafile', 'var')
    load(fullfile(datadir,metafile), 'metadata');
  end
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
  fprintf('Loading data from  %s... ', datapath);
  load(datapath, 'X');

  if isfield(jdat,'filters')
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
  if isfield(jdat,'SanityCheckData')
    disp('PSYCH! This is a simulation.');
    switch jdat.SanityCheckData
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

  if isfield(jdat, 'cvfile')
    cvpath = fullfile(cvdir,cvfile);
    load(cvpath, 'CV');
    outlying_words = reduxFilter.words(filter);
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
    [n,~] = size(X);
    cvind = mod(1:n, ncv)+1;
    cvind = cvind(randperm(n));
  end

  %% ----------------Visual, Audio or Semantic similarities and processing----------------

  if isfield(jdat,'SanityCheckModel')
    switch jdat.SanityCheckModel
    case 'shuffle'
      disp('Shuffling Similarity Matrix!!!')
      simpath = fullfile(datadir,simfile);
      allSimStructs = load(simpath);
      S = allSimStructs.(jdat.simtype);
      shidx = randperm(size(S,1));
      S = S(shidx,shidx);

    case 'random'
      fprintf('Generating totally random Similarity Matrix!!! ')
      simpath = fullfile(datadir,simfile);
      allSimStructs = load(simpath);
      S = allSimStructs.(jdat.simtype);
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
    end
  end
  S = S(filter,filter); clear allSimStructs;
  S = S(~finalholdout, ~finalholdout);

  fprintf('Data loaded and processed.\n');

  %% ---------------------Setting algorithm parameters-------------------------
  [Uz, Sz, nz_rows, p1] = learn_similarity_encoding(S, X, lambda_in, cvind, cvholdout, normalize, Gtype, DEBUG, opts); %#ok<ASGLU>

  fprintf('Saving stuff.....\n');

  save(matfilename,'p1','nz_rows','Sz','Uz');

  fprintf('Done!\n');
end
