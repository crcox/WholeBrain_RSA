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
  nchar = 0;
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
%       fprintf(repmat('\b', 1, nchar));
%       nchar = fprintf('%d',i);
      tmp = RESULTS(i);
      sind = sscanf(PARAMS(i).data,'s%d');
      sidx = find([metadata.subject]==sind);
      % A filter to remove Outliers
      z  = metadata(sidx).rowfilter(1:148);% major hack!
      % other filters already have outliers removed, so we need to put them
      % back to full length.
      finalfilter = false(size(z));
      test = false(size(z));
      finalfilter(z) = tmp.finalfilter; % In fact, this should not be necessary, but it's ``wrong'' in the main code.
      test(z) = tmp.cvfilter;
      train = ~test & ~finalfilter & z;
      
%       train = ~tmp.cvfilter & ~tmp.finalfilter;
      [Ctmp,r] = sqrt_truncate_r(S(~finalfilter,~finalfilter), t);
      C = nan(length(z), r);
      C(~finalfilter,:) = Ctmp; clear Ctmp;
      Cz = nan(length(z), r);
%       Cz = tmp.Cz;
      Cz(z,:) = tmp.Cz; % This filtering should not be necessary
      FroErr1 = norm( C(test,:)-Cz(test,:) ,'fro') / norm( C(test,:) , 'fro');
      FroErr2 = norm( C(train,:)-Cz(train,:) ,'fro') / norm( C(train,:) , 'fro');
      tmp.FroErr1 = FroErr1;
      tmp.FroErr2 = FroErr2;
      R(i) = tmp;
    end
    fprintf('\n');
  end
end
