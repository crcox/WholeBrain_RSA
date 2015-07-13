function R = pairWithCoordinates(results, params, varargin)
  p = inputParser();
  addRequired(p,'results',@isstruct);
  addRequired(p,'params',@isstruct);
  addParameter(p,'DataDir','~/Manchester/WholeBrain_RSA/data/rep',@ischar)
  addParameter(p,'MetaFile','metadata_rep2.mat',@ischar)
  addParameter(p,'CoordFile','~/Manchester/data/mat/fromRick/coords/mni.mat',@ischar)
  parse(p,results,params,varargin{:});

  RESULTS   = p.Results.results;
  PARAMS    = p.Results.params;
  datadir   = p.Results.DataDir;
  metafile  = fullfile(datadir,p.Results.MetaFile);
  coordfile = p.Results.CoordFile;

  load(metafile,'metadata');
  load(coordfile,'xyz');

  DATA = unique({PARAMS.data});
  SID  = cellfun(@(x) sscanf(x, 's%02d'), {PARAMS.data});
  N = length(DATA);
  for i = 1:N
    d = DATA{i};
    sid = SID(i);
    subject = find([metadata.subject] == sid);
    meta = metadata(subject);
    idx = find(strcmp(d, {PARAMS.data}));
    z = meta.colfilter;
    for ii = idx
      tmp = RESULTS(ii);
      nz_rows = RESULTS(ii).nz_rows;
      tmp.xyz = xyz{subject}(nz_rows,:);
      R(ii) = tmp;
    end
  end
end
