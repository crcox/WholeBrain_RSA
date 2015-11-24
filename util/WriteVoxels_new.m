function WriteVoxels_new(outdir,results,params,varargin)
% This function will write out a single directory worth of data, where a
% directory will include all cross validation instances for all subjects for a
% particular configuration of parameters.
% The file naming convention is subject_finalholdout_cvholdout.CoordLabel.
  p = inputParser();
  addRequired(p, 'outdir');
  addRequired(p, 'results');
  addRequired(p, 'params');
  addParameter(p, 'CoordLabel', 'mni');
  addParameter(p, 'SubjectNumberFMT', 's%d_avg.mat');
  parse(p, outdir, results, params, varargin{:});

  outdir    = p.Results.outdir;
  results   = p.Results.results;
  params    = p.Results.params;
  xyzlab    = p.Results.CoordLabel;
  sjfmt     = p.Results.SubjectNumberFMT;

  n = length(results);
  fprintf('Writing files:\n');
  for i = 1:n
    R  = results(i);
    sj = R.subject;
    fh = R.finalholdout;
    cv = R.cvholdout;

    fname = sprintf('%02d_%02d_%02d.%s', sj, fh, cv, xyzlab);
    fpath = fullfile(outdir, fname);
    fprintf('%s\n', fpath);

    idx = find(strcmp(xyzlab,{R.coords.label}));
    xyz = R.coords(idx).xyz;
    dlmwrite(fpath, xyz, ' ');
  end
  fprintf('Done.\n');
end
