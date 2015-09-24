function WriteVoxels(outdir,select,results,params,varargin)
% This function will write out a single directory worth of data, where a
% directory will include all cross validation instances for all subjects for a
% particular configuration of parameters.
% The file naming convention is subject_finalholdout_cvholdout.CoordLabel.
  p = inputParser();
  addRequired(p, 'outdir');
  addRequired(p, 'select')
  addRequired(p, 'results');
  addRequired(p, 'params');
  addParameter(p, 'CoordLabel', 'mni');
  addParameter(p, 'SubjectNumberFMT', 's%d_rep.mat');
  parse(p, outdir, select, results, params, varargin{:});

  outdir    = p.Results.outdir;
  results   = p.Results.results;
  select    = p.Results.select;
  params    = p.Results.params;
  xyzlab    = p.Results.CoordLabel;
  sjfmt     = p.Results.SubjectNumberFMT;

  assert(all(structfun(@isscalar, select)|structfun(@ischar, select)));

  fields = fieldnames(select);
  n = length(fields);
  z = true(size(params));
  for i = 1:n
    key = fields{i};
    val = select.(key);
    if isscalar(val)
      z = z & ([params.(key)] == val);
    elseif ischar(val)
      z = z & strcmp(val, {params.(key)});
    end
  end

  results = results(z);
  params = params(z);

  n = length(params);
  fprintf('Writing files.\n');
  j = 0;
  for i = 1:n
    P  = params(i);
    R  = results(i);
    dt = P.data;
    fh = P.finalholdout;
    cv = P.cvholdout;
    sj = sscanf(dt, sjfmt);

    fname = sprintf('%02d_%02d_%02d.%s', sj, fh, cv, xyzlab);
    fpath = fullfile(outdir, fname);

    idx = find(strcmp(xyzlab,{R.coords.label}));
    if any(R.nz_rows>1)
      v = R.nz_rows(:);
      xyz = [R.coords(idx).xyz, v(v>0)];
    else
      xyz = R.coords(idx).xyz;
    end

    if size(xyz,1) == 0
      if j == 0
        fprintf('* These files will be empty because there are no voxels.\n');
      end
      j = j + 1;
      fprintf('%s\n', fpath);
    end
    dlmwrite(fpath, xyz, ' ');
  end
  if j > 0
    fprintf('Alert: %d empty files written.\n', j);
  end
  fprintf('Done.\n');
end
