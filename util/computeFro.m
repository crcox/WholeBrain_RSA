function R = computeFro(results,params,varargin)
  p = inputParser();
  addRequired(p,'results',@isstruct);
  addRequired(p,'params',@isstruct);
  addParameter(p,'DataDir','~/Manchester/WholeBrain_RSA/data/rep',@ischar)
  addParameter(p,'MetadataFile','metadata_rep.mat',@ischar)
  addParameter(p,'SimilarityFile','model_visual_rep.mat',@ischar)
  addParameter(p,'SimilarityType','earthmover',@ischar)
  parse(p,results,params,varargin{:});

  RESULTS   = p.Results.results;
  PARAMS    = p.Results.params;
  datadir   = p.Results.DataDir;
  metafile  = fullfile(datadir, p.Results.MetadataFile);
  simfile   = fullfile(datadir, p.Results.SimilarityFile);
  simtype   = p.Results.SimilarityType;

  load(metafile);
  AllS = load(simfile);
  S = AllS.(simtype);

  tau = [PARAMS.tau];
  TAU = unique(tau);
  N = length(TAU);
  for ii = 1:N
    t = TAU(ii);
    disp(t)
    tindex = find(tau == t);
    if (t == 0.1)
      continue;
    end
    for i = tindex
      if isempty(RESULTS(i).Cz)
        continue;
      end
      tmp = RESULTS(i);
      sind = sscanf(PARAMS(i).data,'s%d');
      sidx = find([metadata.subject]==sind);
      z  = metadata(sidx).rowfilter(1:148);% major hack!
      test = tmp.cvfilter;
      train = ~tmp.cvfilter & ~tmp.finalfilter;
      C = sqrt_truncate_r(S(z,z), t);
      Cz = tmp.Cz;
      FroErr1 = norm( C(test,:)-Cz(test,:) ,'fro') / norm( C(test,:) , 'fro');
      FroErr2 = norm( C(train,:)-Cz(train,:) ,'fro') / norm( C(train,:) , 'fro');
      tmp.FroErr1 = FroErr1;
      tmp.FroErr2 = FroErr2;
      R(i) = tmp;
    end
  end
end
