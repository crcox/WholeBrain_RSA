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

  n = length(jobdirs);
  nchar = 0;
  fprintf('Loading job ');
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
    R = load(resultpath);
    nrow = size(R.Sz, 1);
    if params(i).bias
      ncol = size(R.Uz, 1)-1;
      R.nz_rows(end) = [];
    else
      ncol = size(R.Uz, 1);
    end

    % Construct filters
    if isempty(filter_labels)
      rowfilter = true(1, nrow);
      colfilter = true(1, ncol);
    else
      if ~iscell(filter_labels);
        filter_labels = {filter_labels};
      end
      % M.filter points to a structured array of filters.
      % First, force filters to a common orientation.
      for ii = 1:numel(M.filter)
        M.filter(ii).filter = forceRowVec(M.filter(ii).filter);
      end
      % Then select the filters
      z = false(1,numel(M.filter));
      for f = filter_labels;
        z(strcmp(f, {M.filter.label})) = true;
      end
      z = z & strcmp(data_varname, {M.filter.subset});

      filters.row = M.filter(z & [M.filter.dimension]==1);
      filters.col = M.filter(z & [M.filter.dimension]==2);
      if isempty(filters.row)
        rowfilter = true(1, nrow);
      else
        rowfilter = all(cat(1, filters.row.filter),1);
      end
      if isempty(filters.col)
        colfilter = true(1, ncol);
      else
        colfilter = all(cat(1, filters.col.filter),1);
      end
      clear filters;
    end
    cvfilter    = cvind == CV(rowfilter,cvscheme);
    finalfilter = finalind == CV(rowfilter,cvscheme);

    % apply filters
    S = SS(rowfilter,rowfilter);
    R.S = S(~finalfilter,~finalfilter);
    R.S_test = S(cvfilter,cvfilter);
    XYZ = XYZ(colfilter,:);

    % log coordinates of selected voxels
    R.coords.label = coord_varname;
    R.coords.xyz = XYZ(R.nz_rows,:);

    % log metadata
    R.cvind = cvind;
    R.cvfilter = cvfilter;
    R.finalfilter = finalfilter;

    % Drop any variables that should not be held in memory
    if ~isempty(SKIP)
      nskip = length(SKIP);
      for ii = 1:nskip
        R = rmfield(R, SKIP{ii});
      end
    end

    results(i) = R;

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
