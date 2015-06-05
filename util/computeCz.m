function R = computeCz(results,params,varargin)
  p = inputParser();
  addRequired(p,'results',@isstruct);
  addRequired(p,'params',@isstruct);
  addParameter(p,'DataDir','~/MCW/WholeBrain_RSA/data/animals',@ischar)
  addParameter(p,'MetadataFile','metadata.mat',@ischar)
  parse(p,results,params,varargin{:});

  RESULTS   = p.Results.results;
  PARAMS    = p.Results.params;
  datadir   = p.Results.DataDir;
  metafile  = fullfile(datadir,p.Results.MetadataFile);

  filters = PARAMS(1).filters;
  load(metafile);
  DATA = unique({params.data});
  N = length(DATA);
  for ii = 1:N
    d = DATA{ii};
    datafile = fullfile(datadir,d);
    assert(exist(datafile,'file')==2,'%s does not exist. Exiting...',datafile);
    fprintf('Loading %s...\n', datafile);
    load(datafile,'X');
    X = filterData(X, metadata, filters);
    dindex = find(strcmp(d, {params.data}));
    n = length(dindex);
    ii = 0;
    nchar=0;
    for i = dindex
      ii = ii + 1;
      fprintf(repmat('\b',1,nchar'));
      nchar = fprintf('%d%%', round((ii/n)*100));
      jobdir = PARAMS(i).jobdir;
      rfile = fullfile(jobdir,'results.mat');
      if ~exist(rfile, 'file')
        continue
      end
      tmp = RESULTS(i);
      cvind = tmp.cvind;
      finalind = PARAMS(i).finalholdout;
      cvind(cvind>finalind) = cvind(cvind>finalind) - 1;
      load(rfile,'Uz');
      tmp.Cz = X*Uz{cvind};
      R(i) = tmp;
    end
    fprintf('\n');
  end

end

function X = filterData(X,metadata,filters)
  if ~isempty(filters)
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
  allzero = any(X); % Identify columns with data
  [~, reduxFilter] = removeOutliers(X);
  % Note: reduxFilter field names are reversed...
  filter = filter & reduxFilter.words';
  vox = allzero & reduxFilter.voxels;
  X = X(filter,vox);
end

