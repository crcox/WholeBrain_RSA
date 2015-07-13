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
  addParameter(p,'SkipFields',[])
  parse(p,varargin{:});

  resultdir = p.Results.ResultDir;
  datadir   = p.Results.DataDir;
  simfile   = p.Results.SimilarityFile;
  simtype   = p.Results.SimilarityType;
  cvsfile   = p.Results.CVSchemesFile;
  metafile  = p.Results.MetadataFile;
  filters   = p.Results.Filters;
  sortjobs  = p.Results.SortJobs;
  SKIP      = p.Results.SkipFields;
  allfiles  = dir(resultdir);
  alldirs   = allfiles([allfiles.isdir]);
  jobdirs   = SelectJobDirs(alldirs, sortjobs);

  if ~isempty(SKIP)
    SkipStr = SKIP;
    SkipStr{end} = ['and ', SKIP{end}];
    SkipStr = strjoin(SkipStr, ', ');
    fprintf('Fields %s will be skipped.\n', SkipStr);
  end

  cvspath = fullfile(datadir,cvsfile);
  load(cvspath, 'CV');

  simpath  = fullfile(datadir,simfile);
  alltypes = load(simpath);
  SS       = alltypes.(simtype);

  metapath = fullfile(datadir,metafile);
  load(metapath, 'metadata');

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

    sind = sscanf(params(i).data,'s%d');
    sidx = find([metadata.subject]==sind);
    m = length(filters);
    for ii = 1:m
        fname = filters{ii};
        z = metadata(sidx).(fname);
        if ii == 1;
            filter = z;
        else
            filter = filter & z;
        end
    end
    % Select appropriate items from SS
    S = SS(filter,filter);

    cvind       = params(i).cvholdout;
    finalind    = params(i).finalholdout;
    cvscheme    = params(i).cvscheme;
    cvfilter    = cvind == CV(filter,cvscheme); %#ok<NODEF>
    finalfilter = finalind == CV(filter,cvscheme);

    if finalind > 0
        if cvind > finalind
          cv = cvind - 1;
        else
          cv = cvind;
        end
    else
        cv = cvind;
    end

    resultfile      = fullfile(jobdir, 'results.mat');
    if exist(resultfile, 'file')
      tmp             = load(resultfile);
      tmp.p1          = selectcv(tmp.p1,cv);
      tmp.p2          = selectcv(tmp.p2,cv);
      tmp.err1        = selectcv(tmp.err1,cv);
      tmp.err2        = selectcv(tmp.err2,cv);
      tmp.cor1        = selectcv(tmp.cor1,cv);
      tmp.cor2        = selectcv(tmp.cor2,cv);
      tmp.Uz          = selectcv(tmp.Uz,cv);
      tmp.Sz          = selectcv(tmp.Sz,cv);
      tmp.nz_rows     = selectcv(tmp.nz_rows,cv);
      tmp.S           = S(~finalfilter,~finalfilter);
      tmp.S_test      = S(cvfilter,cvfilter);
      tmp.Sz_test     = tmp.Sz(cvfilter(~finalfilter),cvfilter(~finalfilter));
      tmp.cvind       = cvind;
      tmp.cvfilter    = cvfilter;
      tmp.finalfilter = finalfilter;
      if ~isempty(SKIP)
        nskip = length(SKIP);
        for ii = 1:nskip
          tmp = rmfield(tmp, SKIP{ii});
        end
      end
      results(i)  = tmp;
    end
  end
  fprintf('\n')
end

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
