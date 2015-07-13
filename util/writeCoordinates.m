function writeCoordinates(results, params)
  p = inputParser();
  addRequired(p,'results',@isstruct);
  addRequired(p,'params',@isstruct);
  parse(p,results,params);

  RESULTS   = p.Results.results;
  PARAMS    = p.Results.params;

  N = length(PARAMS);
  for i = 1:N
    jobdir = PARAMS(i).jobdir;
    nzfile = fullfile(jobdir, 'xyz_nz.txt');
    nzr = RESULTS(i).nz_rows;
    if all(size(nzr)>1)
      nzr = nzr(:,any(nzr));
      counts = sum(nzr);
    else
      counts = ones(nnz(nzr),1);
    end
    out = [RESULTS(i).xyz, counts];
    dlmwrite(nzfile, out, ' ');
  end
end
