function WriteTable(filename,results,varargin)
  p = inputParser();
  addRequired(p, 'filename');
  addRequired(p, 'results');
  addParameter(p, 'fields',{})
  addParameter(p, 'overwrite', false);
  parse(p, filename, results, varargin{:});

  filename  = p.Results.filename;
  results   = p.Results.results;
  fields    = p.Results.fields;
  OVERWRITE = p.Results.overwrite;

  if OVERWRITE
    fid = fopen(filename, 'w');
  else
    if exist(filename, 'file');
      error('%s already exists. To overwrite, use the set overwrite option to true.', filename);
    else
      fid = fopen(filename, 'w');
    end
  end

  % List all possible fields and their desired format code
  fieldFMT = struct( ...
    'iter'         , '%d'   , ...
    'job'          , '%d'   , ...
    'subject'      , '%d'   , ...
    'finalholdout' , '%d'   , ...
    'cvholdout'    , '%d'   , ...
    'data'         , '%s'   , ...
    'target'       , '%s'   , ...
    'Gtype'        , '%s'   , ...
    'regularization', '%s'   , ...
    'lambda'       , '%.4f' , ...
    'lambda1'      , '%.4f' , ...
    'LambdaSeq'    , '%s'   , ...
    'tau'          , '%.4f' , ...
    'normalize'    , '%s'   , ...
    'bias'         , '%d'   , ...
    'RandomSeed'   , '%d'   , ...
    'nVoxel'       , '%d'   , ...
    'p1'           , '%.4f' , ...
    'p2'           , '%.4f' , ...
    'cor1'         , '%.4f' , ...
    'cor2'         , '%.4f' , ...
    'err1'         , '%.4f' , ...
    'err2'         , '%.4f' , ...
    'FroErr1'      , '%.4f' , ...
    'FroErr2'      , '%.4f' , ...
    'nz_rows'      , '%d'   , ...
    'nvox'         , '%d'   , ...
    'nzv'          , '%d');

  hdrFMT = strjoin(repmat({'%s'},1,length(fields)),',');
  tmp = cellfun(@(x) fieldFMT.(x), fields, 'unif', 0);
  dataFMT = strjoin(tmp,',');
  fprintf(fid,[hdrFMT,'\n'],fields{:});

  N = length(results);
  for i = 1:N
    R = results(i);
    out = cell(1,length(fields));
    for j = 1:length(fields);
      key = fields{j};
      if isfield(R,key)
        out{j} = R.(key);
      end
    end
    fprintf(fid,[dataFMT,'\n'], out{:});

  end
  fclose(fid);
end
