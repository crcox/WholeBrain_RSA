function [Params, Results] = HTCondorLoad(ResultDir, varargin)
% HTCondorLoad Load results from a collection of HTCondor jobs.
%   [params] = HTCondorLoad(ResultDir); avoids loading any data.
%   [params, results] = HTCondorLoad(ResultDir);
%   [params, results] = HTCondorLoad(ResultDir, varargin);
%
%   OPTIONS:
%   ResultFile : The name of the file written by the HTCondor job. Will
%                contain a structure.
%   ParamFile  : The name of the json file that defines the job.
%   SortJobs   : Will ensure that the internal list of job directories is
%                sorted by directory name before loading any data.
%   SkipFields : A cell array of field names in the results structure. Any
%                field name listed here will not be loaded. Useful if there
%                are large data elements that you do not need. Controls
%                memory footprint and speeds access.
%   JobList    : A list of job directory names to load.
%
%  DEPENDENCIES:
%  - JSONlab http://www.mathworks.com/matlabcentral/fileexchange/33381
  p = inputParser();
  addRequired(p,'ResultDir', @ischar);
  addParameter(p,'ResultFile','results.mat',@ischar);
  addParameter(p,'ParamFile','params.json',@ischar);
  addParameter(p,'SortJobs',false,@islogical)
  addParameter(p,'SkipFields',[])
  addParameter(p,'JobList',{},@iscellstr)
  addParameter(p,'legacy',false,@islogical)
  addParameter(p,'quiet',false,@islogical)
  parse(p,ResultDir, varargin{:});

  RESULT_DIR  = p.Results.ResultDir;
  RESULT_FILE = p.Results.ResultFile;
  PARAMS_FILE = p.Results.ParamFile;
  SORT_JOBS   = p.Results.SortJobs;
  SKIP        = p.Results.SkipFields;
  QUIET       = p.Results.quiet;
  LEGACY      = p.Results.legacy;
  jobDirs     = p.Results.JobList;

  if isempty(jobDirs)
      allFiles = dir(RESULT_DIR);
      allDirs  = allFiles([allFiles.isdir]);
      jobDirs  = SelectJobDirs(allDirs, PARAMS_FILE, RESULT_DIR, SORT_JOBS);
  end
  nJobDirs = numel(jobDirs);

  if ~isempty(SKIP)
      SkipStr = SKIP;
      SkipStr{end} = ['and ', SKIP{end}];
      SkipStr = strjoin(SkipStr, ', ');
      fprintf('Fields %s will be skipped.\n', SkipStr);
  end

%% Preallocate result structure
% Assumption is that all result files are identical size with same fields
  i = 0;
  fileExists = false;
  if nargout > 1
    while ~fileExists
        i = i + 1;
        jobDir     = fullfile(RESULT_DIR, jobDirs{i});
        resultPath = fullfile(jobDir, RESULT_FILE);
        fileExists = exist(resultPath, 'file') == 2;
    end
    paramsPath = fullfile(jobDir, PARAMS_FILE);
    tmp_p = loadjson(paramsPath);
    tmp_p.jobdir = 0;
    tmp_p.subject = 0;
    Params(nJobDirs) = tmp_p;
    tmp = load(resultPath);
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
    N = nJobDirs * numel(R);
    Results = R(1);
    Results.job = 0;
    Results(N).job = 0;
  end

  nchar = 0;
  if ~QUIET
      fprintf('Loading job ');
  end
  cursor = 0;
  for i = 1:nJobDirs;
    if ~QUIET
        fprintf(repmat('\b', 1, nchar));
        nchar = fprintf('%d of %d', i, nJobDirs);
    end

    % load parameter file
    jobDir      = fullfile(RESULT_DIR, jobDirs{i});
    paramsPath  = fullfile(jobDir,PARAMS_FILE);
    tmp         = loadjson(paramsPath);
    tmp.jobdir  = jobDir;
    tmp.subject = sscanf(tmp.data,'s%02d');
    Params(i)   = tmp;
    clear tmp;

    if nargout > 1
      % load results file
      resultPath = fullfile(jobDir, RESULT_FILE);
      if ~exist(resultPath, 'file')
          continue;
      end
      tmp = load(resultPath);
      if isfield(tmp,'results')
        R = tmp.results;
      else
        R = tmp;
      end
      if LEGACY
        tmp = [fieldnames(R)';cellfun(@(x) mat2cell(x, ones(size(x,1),1), size(x,2)), struct2cell(R), 'Unif', 0)'];
        R = struct(tmp{:});
      end
      [R.job] = deal(i);
      if ~isempty(SKIP)
          R = rmfield(R, SKIP);
      end
      a = cursor + 1;
      b = cursor + numel(R);
      if LEGACY
        Params(i).cvholdout = mat2cell(Params(i).cvholdout,1,ones(1,numel(Params(i).cvholdout)));
        Params(i).filters = {Params(i).filters};
        Params(i).COPY = {Params(i).COPY};
        tmp = [fieldnames(Params)';struct2cell(Params(i))'];
        ParamsAB = struct(tmp{:});

        fnames = fieldnames(Results);
        nf = numel(fnames);
        for iField = 1:nf
          fn = fnames{iField};
          iiR = 0;
          for iR = a:b
            iiR = iiR + 1;
            if isfield(R,fn)
              Results(iR).(fn) = R(iiR).(fn);
            elseif isfield(Params,fn)
              Results(iR).(fn) = ParamsAB(iiR).(fn);
            else
              % do nothing
            end
          end
        end
      else
        Results(a:b) = R;
      end
      cursor = b;
      a = cursor + 1;
      if a < b
          Results(a:b) = [];
      end
    end
  end
  if ~QUIET
      fprintf('\n')
  end
end

function jobDirs = SelectJobDirs(dirs, paramsFile, resultDir, sort)
  p = inputParser;
  addRequired(p, 'dirs');
  addRequired(p, 'ParamsFile', @ischar);
  addRequired(p, 'ResultDir', @ischar);
  addRequired(p, 'sort', @islogical);
  parse(p, dirs, paramsFile, resultDir, sort);

  JOB_DIRS = p.Results.dirs;
  PARAMS_FILE = p.Results.ParamsFile;
  RESULT_DIR = p.Results.ResultDir;
  SORT = p.Results.sort;

  N = length(JOB_DIRS);
  isJobDir = false(N,1);
  for ii = 1:N
    jobDir = fullfile(RESULT_DIR,JOB_DIRS(ii).name);
    % Check if a special dir (current or parent)
    if any(strcmp(jobDir,{'.','..'}))
      continue
    end
    % Check if contains parameter file.
    paramsPath = fullfile(jobDir, PARAMS_FILE);
    if exist(paramsPath, 'file')
      isJobDir(ii) = true;
    end
  end
  jobDirs = JOB_DIRS(isJobDir);

  if SORT
    try
      jobs = cellfun(@(x) sscanf(x,'%d'), {jobDirs.name});
      disp('Sorting jobs numerically.')
    catch
      jobs = {jobDirs.name};
      disp('Sorting jobs alphabetically.')
    end
    [~,ix] = sort(jobs);
    jobDirs = jobDirs(ix);
  end
  jobDirs = {jobDirs.name};
  fprintf('Found %d job directories.\n', length(jobDirs))
end
