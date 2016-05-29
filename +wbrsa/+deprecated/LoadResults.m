function [results, params] = LoadResults(varargin)
  p = inputParser();
  addParameter(p,'ResultDir','.',@ischar);
  addParameter(p,'DataDir','~/data/Manchester/WholeBrain_RSA/data/avg',@ischar)
  addParameter(p,'MetadataFile','metadata_avg.mat',@ischar)
  addParameter(p,'MetadataVarname','metadata',@ischar)
  addParameter(p,'ResultFile','results.mat',@ischar);
  addParameter(p,'ParamFile','params.json',@ischar);
  addParameter(p,'SortJobs',false,@islogical)
  addParameter(p,'SkipFields',[])
  parse(p,varargin{:});

  resultdir     = p.Results.ResultDir;
  datadir       = p.Results.DataDir;
  metafile      = p.Results.MetadataFile;
  meta_varname  = p.Results.MetadataVarname;
  sortjobs      = p.Results.SortJobs;
  resultfile    = p.Results.ResultFile;
  paramsfilename = p.Results.ParamFile;
  SKIP          = p.Results.SkipFields;

  allfiles      = dir(resultdir);
  alldirs       = allfiles([allfiles.isdir]);
  jobdirs       = SelectJobDirs(alldirs, 'ParamsFilename', paramsfilename, 'root',resultdir, 'sort', sortjobs);

  if ~isempty(SKIP)
    SkipStr = SKIP;
    SkipStr{end} = ['and ', SKIP{end}];
    SkipStr = strjoin(SkipStr, ', ');
    fprintf('Fields %s will be skipped.\n', SkipStr);
  end

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
    paramspath  = fullfile(jobdir,paramsfilename);
    tmp         = loadjson(paramspath);
    tmp.jobdir  = jobdir;
    params(i)   = tmp;
    clear tmp;

    % load results file
    resultpath = fullfile(jobdir, resultfile);
    if ~exist(resultpath, 'file')
      continue;
    end
    tmp = load(resultpath);
    % For back compatibility
    if ~isfield(tmp,'results')
      tmp = repackageresults(tmp);
    end
    R = tmp.results;
    [R.job] = deal(i);
    if ~isempty(SKIP)
      R = rmfield(R, SKIP);
    end
    if ~isfield(results, 'data_varname') && isfield(R,'data_varname')
      [results.data_varname] = deal([]);
    end
    if ~isfield(results, 'filters')
      [results.filters] = deal([]);
    end
    if isfield(params, 'data_varname') && ~isfield(R, 'data_varname')
      [R.data_varname] = deal(params(i).data_varname);
    end
    if isfield(params, 'filters') && ~isfield(R, 'filters')
      [R.filters] = deal(params(i).filters);
    end
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
function R = repackageresults(r)
  m = numel(r.err1);
  R.results = init_results();
  R.results(m).job = [];
  nz_rows = mat2cell(r.nz_rows, ones(m,1), size(r.nz_rows,2));
  p1 = mat2cell(r.p1, ones(m,1), 1);
  p2 = mat2cell(r.p2, ones(m,1), 1);
  cor1 = mat2cell(r.cor1, ones(m,1), 1);
  cor2 = mat2cell(r.cor2, ones(m,1), 1);
  p1t = mat2cell(r.p1t, ones(m,1), 1);
  p2t = mat2cell(r.p2t, ones(m,1), 1);
  cor1t = mat2cell(r.cor1t, ones(m,1), 1);
  cor2t = mat2cell(r.cor2t, ones(m,1), 1);
  err1 = mat2cell(r.err1, ones(m,1), 1);
  err2 = mat2cell(r.err2, ones(m,1), 1);
  iter = r.iter;
  [R.results.nz_rows] = deal(nz_rows{:});
  [R.results.p1] = deal(p1{:});
  [R.results.p2] = deal(p2{:});
  [R.results.cor1] = deal(cor1{:});
  [R.results.cor2] = deal(cor2{:});
  [R.results.p1t] = deal(p1t{:});
  [R.results.p2t] = deal(p2t{:});
  [R.results.cor1t] = deal(cor1t{:});
  [R.results.cor2t] = deal(cor2t{:});
  [R.results.err1] = deal(err1{:});
  [R.results.err2] = deal(err2{:});
  [R.results.iter] = deal(iter);
end

function jobdirs = SelectJobDirs(dirs,varargin)
  p = inputParser;
  addRequired(p, 'dirs');
  addParameter(p, 'ParamsFilename','params.json', @ischar);
  addParameter(p, 'root','.', @ischar);
  addParameter(p, 'sort',false, @islogical);
  parse(p, dirs, varargin{:});

  dirs = p.Results.dirs;
  paramsfilename = p.Results.ParamsFilename;
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
    paramspath = fullfile(jobdir, paramsfilename);
    if exist(paramspath, 'file')
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
