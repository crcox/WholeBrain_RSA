function RESULTS = PairWithCoordinates(results, metadata, coords, coordslabel, varargin)
  p = inputParser();
  addRequired(p,'results',@isstruct);
  addRequired(p,'metadata',@isstruct)
  addRequired(p,'coords',@isstruct)
  addRequired(p,'coordslabel',@ischar)
  parse(p,results, metadata, coords, coordslabel, varargin{:});

  RESULTS  = p.Results.results;
  metadata = p.Results.metadata;
  coords = p.Results.coords;
  coordslabel  = p.Results.coordslabel;

  ss = [results.subject];
  N = length(ss);
  for i = 1:N
    subject = ss(i);
    idx = find([results.subject] == subject);
    BiasUnit = results(idx(1)).bias;
    z = [coords.subject] == subject;
    if nnz(z) > 1
      error('Subject numbers are not unique in coords struct.\n');
    end
    xyz0 = coords(z).(coordslabel);
    for ii = idx
      tmp = RESULTS(ii);
      if BiasUnit
        nz_rows = RESULTS(ii).nz_rows(1:end-1);
      else
        nz_rows = RESULTS(ii).nz_rows;
      end
      dvn = results(ii).data_varname;
      czlabs = results(ii).filters;
      colfilter = true(size(xyz0, 1), 1);
      if iscell(czlabs)
        for i = 1:numel(czlabs)
          lab = czlabs{i};
          z = strcmp(dvn, {metadata(subject).filter.subset}) & ...
              strcmp(lab, {metadata(subject).filter.label}) & ...
              [metadata(subject).filter.dimension] == 2;
          if any(z)
            f = metadata(subject).filter(z).filter;
            colfilter = colfilter(:) & f(:);
          end
        end
      end
      xyz = xyz0(colfilter, :);
      RESULTS(ii).coords.label = coordslabel;
      RESULTS(ii).coords.xyz = xyz(nz_rows,:);
    end
  end
end
