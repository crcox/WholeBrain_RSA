function [results, params] = LoadResults(varargin)
  p = inputParser();
  addParameter(p,'ResultDir','.',@ischar);
  addParameter(p,'DataDir','~/MCW/WholeBrain_RSA/data/animals',@ischar)
  addParameter(p,'SimilarityFile','semantic_model.mat',@ischar)
  addParameter(p,'SimilarityType','cor',@ischar)
  addParameter(p,'CVSchemesFile','CV_schemes.mat',@ischar)
  addParameter(p,'MetadataFile','metadata.mat',@ischar)
  addParameter(p,'Filters',{'TrueAnimals'},@iscell)
  addParameter(p,'SortJobs',false,@islogical)
  parse(p,varargin{:});

  resultdir = p.Results.ResultDir;
  datadir   = p.Results.DataDir;
  simfile   = p.Results.SimilarityFile;
  simtype   = p.Results.SimilarityType;
  cvsfile   = p.Results.CVSchemesFile;
  metafile  = p.Results.MetadataFile;
  filters   = p.Results.Filters;
  sortjobs  = p.Results.SortJobs;
  allfiles  = dir(resultdir);
  alldirs   = allfiles([allfiles.isdir]);
  jobdirs   = SelectJobDirs(alldirs, sortjobs);

  cvspath = fullfile(datadir,cvsfile);
  load(cvspath, 'CV');

  simpath  = fullfile(datadir,simfile);
  alltypes = load(simpath);
  SS       = alltypes.(simtype);

  metapath = fullfile(datadir,metafile);
  load(metapath, 'metadata');
  n = length(filters);
  for i = 1:n
      sname = filters{i};
      z = metadata.(sname);
      if i == 1;
          filter = z;
      else
          filter = filter & z;
      end
  end

  % Select appropriate items from SS
  S = SS(filter,filter);

  n = length(jobdirs);
  nchar = 0;
  fprintf('Loading job ');
  for i = 1:n;
    fprintf(repmat('\b', 1, nchar));
    nchar = fprintf('%d of %d', i, n);

    jobdir      = jobdirs(i).name;
    paramfile   = fullfile(jobdir,'params.json');
    tmp         = loadjson(paramfile);
    tmp.jobdir  = jobdir;
    params(i)   = tmp;

    cvind       = params(i).cvholdout;
    finalind    = params(i).finalholdout;
    cvscheme    = params(i).cvscheme;
    cvfilter    = cvind == CV(filter,cvscheme); %#ok<NODEF>
    finalfilter = finalind == CV(filter,cvscheme);

    if cvind > finalind
      cv = cvind - 1;
    else
      cv = cvind;
    end

    resultfile      = fullfile(jobdir, 'results.mat');
    tmp             = load(resultfile);
    tmp.p1          = tmp.p1(cv);
    tmp.p2          = tmp.p2(cv);
    tmp.err1        = tmp.err1(cv);
    tmp.err2        = tmp.err2(cv);
    tmp.cor1        = tmp.cor1(cv);
    tmp.cor2        = tmp.cor2(cv);
    tmp.Uz          = tmp.Uz{cv};
    tmp.Sz          = tmp.Sz{cv};
    tmp.nz_rows     = tmp.nz_rows(cv,:);
    tmp.S           = S(~finalfilter,~finalfilter);
    tmp.S_test      = S(cvfilter,cvfilter);
    tmp.Sz_test     = tmp.Sz(cvfilter(~finalfilter),cvfilter(~finalfilter));
    tmp.cvind       = cvind;
    tmp.cvfilter    = cvfilter;
    tmp.finalfilter = finalfilter;
    results(i)  = tmp;
  end
  fprintf('\n')
end
function jobdirs = SelectJobDirs(dirs,SORT)
  N = length(dirs);
  isJobDir = false(N,1);
  for ii = 1:N
    jobdir = dirs(ii).name;
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
