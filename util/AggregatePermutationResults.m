function PR = AggregatePermutationResults(results, varargin)
  p = inputParser();
  addRequired(p,'results', @isstruct);
  parse(p, results, varargin{:});

  results = p.Results.results;

  ss = [results.subject];
  subject = unique(ss);
  N = numel(subject);
  PR(N) = init_results('Template','permtest');
  permfields = fieldnames(PR);
  for i = 1:N
    s = subject(i);
    idx = find(ss==s);
    BiasUnit = results(idx(1)).bias;

    % Count nz rows in Uz
    nzr = [results(idx).nz_rows];
    if BiasUnit
      nnzr = sum(nzr(1:(end-1),:),2);
    else
      nnzr = sum(nzr(1:end,:),2);
    end

    % Compute and aggregate node strengths
    nodestrength = zeros(1, numel(nnzr));
    nzv = 0;
    for ii = 1:numel(idx);
      j = idx(ii);
      if BiasUnit
        Uz = results(j).Uz(1:(end-1),:);
        z = results(j).nz_rows(1:(end-1),:);
      else
        Uz = results(j).Uz;
        z = results(j).nz_rows;
      end
      Uz = Uz(z,:);
      nzv = nzv + nnz(z);
      nodestrength(z) = nodestrength(z) + sum(abs(Uz * Uz'));
    end
    nodestrength = nodestrength / numel(idx);
    nzv = nzv / numel(idx);

    for ii = 1:numel(permfields);
      f = permfields{ii};
      if isfield(results, f)
        PR(i).(f) = results(idx(1)).(f);
      end
    end
    PR(i).count_nz_rows = nnzr;
    if BiasUnit
      PR(i).nz_rows = [nnzr>0;true];
    else
      PR(i).nz_rows = nnzr>0;
    end
    PR(i).mean_nodestrength = nodestrength;
    PR(i).mean_nzv = nzv;
  end
end
