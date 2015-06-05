function R = computeFro(results,params,varargin)
  p = inputParser();
  addRequired(p,'results',@isstruct);
  addRequired(p,'params',@isstruct);
  addParameter(p,'DataDir','~/MCW/WholeBrain_RSA/data/animals',@ischar)
  addParameter(p,'SimilarityFile','semantic_model.mat',@ischar)
  addParameter(p,'SimilarityType','cor',@ischar)
  parse(p,results,params,varargin{:});

  RESULTS   = p.Results.results;
  PARAMS    = p.Results.params;
  datadir   = p.Results.DataDir;
  simfile   = fullfile(datadir, p.Results.SimilarityFile);
  simtype   = p.Results.SimilarityType;

  AllS = load(simfile);
  S = AllS.(simtype);

  tau = [PARAMS.tau];
  TAU = unique(tau);
  N = length(TAU);
  for ii = 1:N
    t = TAU(ii);
    disp(t)
    C = sqrt_truncate_r(S, t);
    tindex = find(tau == t);
    for i = tindex
      if isempty(RESULTS(i).Cz)
        continue;
      end
      Cz = RESULTS(i).Cz;
      whos C Cz
      tmp = RESULTS(i);
      FroErr = norm( C-Cz ,'fro') / norm( C , 'fro');
      tmp.FroErr = FroErr;
      R(i) = tmp;
    end
  end
end
