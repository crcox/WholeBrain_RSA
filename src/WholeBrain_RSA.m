function WholeBrain_RSA()
  %@AVS = 'auditory', 'visual', 'semantic'
  %@lambda_in = scalar for group lasso and vector for group OWL (inversely proportional to sparsity)
  %@normalize = 1 for normalized feature vectors
  %@subject = '03' for subject 03

  % ----------------------Set parameters-----------------------------------------------
  jdat = loadjson('params.json')
  DEBUG = jdat.debug;
  Gtype = jdat.Gtype;
  lambda_in = jdat.lambda;
  normalize = jdat.normalize
  AVS = jdat.reptype;
  targets = jdat.targets;
  datafile = jdat.data;
  cvfile = jdat.cvfile;
  metafile = jdat.metadata;
  if isfield(jdat,'AdlasOpts')
    opts = jdat.AdlasOpts;
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
      [ddir,dname,ext] = fileparts(datafile);
      datadir = fullfile(root,ddir);
      datafile = strcat(dname,ext);
      [mdir,mname,ext] = fileparts(metafile);
      datadir = fullfile(root,mdir);
      metafile = strcat(mname,ext);
      [cdir,cvname,ext] = fileparts(cvfile);
      cvdir = fullfile(root,cdir);
      cvfile = strcat(cvname,ext);
      assert(strcmp(mdir,ddir) && strcmp(cdir,mdir));
    end

  case 'chris'
    root = '/home/chris/MCW/WholeBrain_RSA';
    datadir = fullfile(root,'data');

  case 'urvashi'
    root = '/mnt/ws/home/uoswal/words/';
    datadir = fullfile(root,'data');

  otherwise
    error('Environment %s not implemented.', jdat.environment);

  end

  load(fullfile(datadir,metafile), 'metadata');
  [path,fname,ext] = fileparts(datafile);
  subjid = sscanf(fname, 's%d');
  subject = find(subjid == [metadata.subject]);

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
    ind = metadata(subject).TrueAnimals | metadata(subject).TrueArtifacts;

  otherwise
    error('Invalid target set: %s.', targets);

  end
  datapath = fullfile(datadir,datafile);
  fprintf('Loading data from  %s... ', datapath);

  % load input files and filter outliers more than 5 std dev away
  load(datapath);
  allzero = any(X); % Identify columns with data
  [~, reduxFilter] = removeOutliers(X);
  % Note: reduxFilter field names are reversed...
  ind = ind & reduxFilter.words';
  vox = allzero & reduxFilter.voxels;
  X = X(ind,vox);

  if isfield(jdat, 'cvfile')
    cvpath = fullfile(cvdir,cvfile);
    load(cvpath, 'CV');
    outlying_words = reduxFilter.words(ind);
    cvind = CV(outlying_words, jdat.cvscheme);
    holdout = cvind == jdat.cvholdout;
    X(holdout,:) = [];
    cvind(holdout) = [];
    cvind(cvind>jdat.cvholdout) = cvind(cvind>jdat.cvholdout) - 1;

  else
    ncv = jdat.ncv;
    [n,d] = size(X);
    cvind = mod(1:n, ncv)+1;
    cvind = cvind(randperm(n));
  end
  fprintf('done.\n');

  %% ----------------Visual, Audio or Semantic similarities and processing----------------

  switch AVS
  case 'auditory'
    error('AVS %s: Not implemented.', AVS)

  case 'visual'
    error('AVS %s: Not implemented.', AVS)

  case 'semantic'
    load(fullfile(datadir,'semantic_model.mat'),'inner'); % inner product of word2vec data
    sem = inner; clear inner;
    S = sem(ind,ind); clear sem;
    S = S(~holdout, ~holdout);

  otherwise
    error('AVS %s not implemented. Should be aud, vis or sem.\n', AVS)

  end

  fprintf('Data loaded and processed.');

  %% ---------------------Setting algorithm parameters-------------------------
  [Uz, nz_rows, p1] = learn_similarity_encoding(S, X, lambda_in, cvind, normalize, Gtype, DEBUG, opts);

  fprintf('Saving stuff.....\n');

  if normalize == 1
    summaryfile = fullfile(root,sprintf('%s_norm_%s_%s.txt',AVS,Gtype));
  else
    summaryfile = fullfile(root,sprintf('%s_%s_%s.txt',AVS,Gtype));
  end

  save(matfilename,'p1','nz_rows','Uz');

  %fid = fopen(summaryfile,'a+');
  %fprintf(fid,'%s: corr: %f %f test error: %f %f with lambda = %f\n',subject,mean(p1),2*std(p1)/sqrt(tune), mean(test_err),2*std(test_err)/sqrt(tune),lambda_min);
  %fclose(fid);

  %fprintf('Final estimate of mean correlation and confidence for subject %s: %f %f \n',subject,mean(p11),2*std(p11)/sqrt(tune));
  %fprintf('Final estimate of mean error in predicted similarity and confidence for subject %s: %f %f \n',subject,mean(test_err),2*std(test_err)/sqrt(tune));

  fprintf('Done!\n');
end
