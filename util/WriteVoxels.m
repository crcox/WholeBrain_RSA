function WriteVoxels(outdir,results,varargin)
% This function will write out a single directory worth of data, where a
% directory will include all cross validation instances for all subjects for a
% particular configuration of parameters.
% The file naming convention is subject_finalholdout_cvholdout.CoordsLabel.
  p = inputParser();
  addRequired(p, 'outdir', @ischar);
  addRequired(p, 'results', @isstruct);
  addParameter(p, 'CoordsLabel', 'mni', @ischar);
  addParameter(p, 'VoxelValue', '', @ischar);
  parse(p, outdir, results, varargin{:});

  outdir    = p.Results.outdir;
  results   = p.Results.results;
  xyzlab    = p.Results.CoordsLabel;
  valuefield = p.Results.VoxelValue;

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
    if ~isempty(valuefield) && isfield(R, valuefield)
      values = R.(valuefield);
      if isnumeric(values)
        values = values(values~=0);
      end
      if size(xyz,1) == numel(values)
        if islogical(values)
          dlmwrite(fpath, xyz(values,:), ' ');
        else
          dlmwrite(fpath, [xyz,values(:)], ' ');
        end
      end
    else
      dlmwrite(fpath, xyz, ' ');
    end
  end
  fprintf('Done.\n');
end
