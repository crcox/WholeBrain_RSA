function [results, params] = LoadResults(varargin)
  p = inputParser();
  addParameter(p,'ResultDir','.',@ischar);
  addParameter(p,'DataDir','~/MCW/WholeBrain_RSA/data/avg',@ischar)
  addParameter(p,'Coordinates','coords.mat')
  addParameter(p,'coord_varname','mni')
  addParameter(p,'SimilarityFile','semantic_model.mat',@ischar)
  addParameter(p,'SimilarityType','cor',@ischar)
  addParameter(p,'CVSchemesFile','CV_schemes.mat',@ischar)
  addParameter(p,'MetadataFile','metadata.mat',@ischar)
  addParameter(p,'Filters',{'TrueAnimals'},@iscell)
  addParameter(p,'ResultFile','',@ischar);
  addParameter(p,'ParamFile','',@ischar);
  addParameter(p,'SortJobs',false,@islogical)
  addParameter(p,'SkipFields',[])
  parse(p,varargin{:});

  resultdir     = p.Results.ResultDir;
  datadir       = p.Results.DataDir;
  simfile       = p.Results.SimilarityFile;
  sim_varname   = p.Results.SimilarityType;
  coordfile     = p.Results.Coordinates;
  coord_varname = p.Results.coord_varname;;
  cvsfile       = p.Results.CVSchemesFile;
  cv_varname    = 'CV';
  metafile      = p.Results.MetadataFile;
  meta_varname  = 'metadata';
  filter_labels = p.Results.Filters;
  sortjobs      = p.Results.SortJobs;
  resultfile    = p.Results.ResultFile;
  paramfile     = p.Results.ParamFile;
  SKIP          = p.Results.SkipFields;
  allfiles      = dir(resultdir);
  alldirs       = allfiles([allfiles.isdir]);
  jobdirs       = SelectJobDirs(alldirs, 'root',resultdir, 'sort', sortjobs);

  if ~isempty(SKIP)
    SkipStr = SKIP;
    SkipStr{end} = ['and ', SKIP{end}];
    SkipStr = strjoin(SkipStr, ', ');
    fprintf('Fields %s will be skipped.\n', SkipStr);
  end

  cvspath = fullfile(datadir,cvsfile);
  StagingContainer = load(cvspath, cv_varname);
  CV = StagingContainer.(cv_varname);

  coordpath = fullfile(datadir,coordfile);
  StagingContainer = load(coordpath, 'coords');
  coords = StagingContainer.('coords');

  simpath  = fullfile(datadir,simfile);
  StagingContainer = load(simpath, sim_varname);
  SS = StagingContainer.(sim_varname);

  metapath = fullfile(datadir, metafile);
  StagingContainer = load(metapath, meta_varname);
  metadata = StagingContainer.(meta_varname);

  %% Preallocate result structure. Assumption is that all result files are
  %identical size with same fields.
  i = 0;
  FileExists = false;
  while ~FileExists
    i = i + 1;
    jobdir     = fullfile(resultdir, jobdirs(i).name);
    resultpath = fullfile(jobdir, resultfile);
    FileExists = exist(resultpath, 'file') == 2;
  end
  tmp = load(resultpath);
  if ~isfield(tmp,'results')
    tmp2 = tmp; clear tmp;
    m = numel(tmp2.err1);
    tmp.results = init_results();
    tmp.results(m).job = [];
    clear tmp2;
  end
  R = tmp.results;
  if ~isempty(SKIP)
    R = rmfield(R, SKIP);
  end
  N = numel(jobdirs) * numel(R);
  results = R(1);
  results.job = 0;
  results(N).job = 0;

  n = length(jobdirs);
  nchar = 0;
  fprintf('Loading job ');
  cursor = 0;
  for i = 1:n;
    fprintf(repmat('\b', 1, nchar));
    nchar = fprintf('%d of %d', i, n);

    % load parameter file
    jobdir      = fullfile(resultdir, jobdirs(i).name);
    parampath   = fullfile(jobdir,paramfile);
    tmp         = loadjson(parampath);
    tmp.jobdir  = jobdir;
    params(i)   = tmp;
    clear tmp;

    data_varname = params(i).data_varname;
    sind = sscanf(params(i).data,'s%d');
    sidx = find([coords.subject]==sind);
    XYZ = coords(sidx).(coord_varname);
    sidx = find([metadata.subject]==sind);
    M = metadata(sidx);

    % Pull cv indexes
    cvind       = params(i).cvholdout;
    finalind    = params(i).finalholdout;
    cvscheme    = params(i).cvscheme;
    if finalind > 0
      if cvind > finalind
        cv = cvind - 1;
      else
        cv = cvind;
      end
    else
        cv = cvind;
    end

    % load results file
    resultpath = fullfile(jobdir, resultfile);
    if ~exist(resultpath, 'file')
      continue;
    end
    tmp = load(resultpath);
    % For back compatibility
    if ~isfield(tmp,'results')
      tmp2 = tmp; clear tmp;
      m = numel(tmp2.err1);
      tmp.results = init_results();
      tmp.results(m).job = [];
      nz_rows = mat2cell(tmp2.nz_rows, ones(m,1), size(tmp2.nz_rows,2));
      p1 = mat2cell(tmp2.p1, ones(m,1), 1);
      p2 = mat2cell(tmp2.p2, ones(m,1), 1);
      cor1 = mat2cell(tmp2.cor1, ones(m,1), 1);
      cor2 = mat2cell(tmp2.cor2, ones(m,1), 1);
      p1t = mat2cell(tmp2.p1t, ones(m,1), 1);
      p2t = mat2cell(tmp2.p2t, ones(m,1), 1);
      cor1t = mat2cell(tmp2.cor1t, ones(m,1), 1);
      cor2t = mat2cell(tmp2.cor2t, ones(m,1), 1);
      err1 = mat2cell(tmp2.err1, ones(m,1), 1);
      err2 = mat2cell(tmp2.err2, ones(m,1), 1);
      iter = tmp2.iter;
      [tmp.results.nz_rows] = deal(nz_rows{:});
      [tmp.results.p1] = deal(p1{:});
      [tmp.results.p2] = deal(p2{:});
      [tmp.results.cor1] = deal(cor1{:});
      [tmp.results.cor2] = deal(cor2{:});
      [tmp.results.p1t] = deal(p1t{:});
      [tmp.results.p2t] = deal(p2t{:});
      [tmp.results.cor1t] = deal(cor1t{:});
      [tmp.results.cor2t] = deal(cor2t{:});
      [tmp.results.err1] = deal(err1{:});
      [tmp.results.err2] = deal(err2{:});
      [tmp.results.iter] = deal(iter);
      clear tmp2;
    end
    R = tmp.results;
    if ~isempty(SKIP)
      R = rmfield(R, SKIP);
    end
    [R.job] = deal(i);
    a = cursor + 1;
    b = cursor + numel(R);
    results(a:b) = R;
    cursor = b;
  end
  a = cursor + 1;
  n = N;
  if a < b
    results(a:b) = [];
  end
  fprintf('\n')
end

%% Private Functions
function y = selectcv(x,cv)
  if numel(x)==1
    if iscell(x);
      y = x{1};
    else
      y = x;
    end
    return
  end

  dim = size(x);
  if iscell(x)
    if length(dim)>1 && all(dim>1)
      y = x(cv,:);
    else
      y = x{cv};
    end
  else
    if length(dim)>1 && all(dim>1)
      y = x(cv,:);
    else
      y = x(cv);
    end
  end
end

function jobdirs = SelectJobDirs(dirs,varargin)
  p = inputParser;
  addRequired(p, 'dirs');
  addParameter(p, 'root','.', @ischar);
  addParameter(p, 'sort',false, @islogical);
  parse(p, dirs, varargin{:});
  
  dirs = p.Results.dirs;
  SORT = p.Results.sort;
  root = p.Results.root;
  
  N = length(dirs);
  isJobDir = false(N,1);
  for ii = 1:N
    jobdir = fullfile(root,dirs(ii).name);
    % Check if a special dir (current or parent)
    if any(strcmp(jobdir,{'.','..'}))
      continue
    end
    % Check if contains parameter file.
    paramfile = fullfile(jobdir, 'params.json');
    if exist(paramfile, 'file')
      isJobDir(ii) = true;
    end
  end
  jobdirs = dirs(isJobDir);

  if SORT
    try
      jobs = cellfun(@(x) sscanf(x,'%d'), {jobdirs.name});
      disp('Sorting jobs numerically.')
    catch
      jobs = {jobdirs.name};
      disp('Sorting jobs alphabetically.')
    end
    [~,ix] = sort(jobs);
    jobdirs = jobdirs(ix);
  end
  fprintf('Found %d job directories.\n', length(jobdirs))
end

function r = forceRowVec(v)
  r = v(1:end);
end

function c = forceColVec(v)
  c = v(:);
end

function results = init_results()
  results.Uz = [];
  results.Cz = [];
  results.Sz = [];
  results.nz_rows = [];
  results.subject = [];
  results.cvholdout = [];
  results.finalholdout = [];
  results.lambda = [];
  results.lambda1 = [];
  results.LambdaSeq = [];
  results.Gtype = [];
  results.bias = [];
  results.normalize = [];
  results.nzv = [];
  results.p1 = [];
  results.p2 = [];
  results.cor1 = [];
  results.cor2 = [];
  results.p1t = [];
  results.p2t = [];
  results.cor1t = [];
  results.cor2t = [];
  results.err1 = [];
  results.err2 = [];
  results.iter = [];
  results.job = [];
end
