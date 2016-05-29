function results = init_results(varargin)
  p = inputParser();
  addParameter(p,'Template','default',@ischar)
  parse(p, varargin{:});

  switch p.Results.Template
  case 'default'
    results.Uz = [];
    results.Cz = [];
    results.Sz = [];
    results.nz_rows = [];
    results.subject = [];
    results.cvholdout = [];
    results.finalholdout = [];
    results.lambda = [];
    results.lambda1 = [];
    results.LambdaSeq = [];
    results.Gtype = [];
    results.bias = [];
    results.normalize = [];
    results.data_varname = [];
    results.filters = [];
    results.nzv = [];
    results.p1 = [];
    results.p2 = [];
    results.cor1 = [];
    results.cor2 = [];
    results.p1t = [];
    results.p2t = [];
    results.cor1t = [];
    results.cor2t = [];
    results.err1 = [];
    results.err2 = [];
    results.iter = [];
    results.job = [];

  case 'permtest'
    results.nodestrength = [];
    results.mean_nodestrength = [];
    results.count_nz_rows = [];
    results.nz_rows= [];
    results.nzv = [];
    results.mean_nzv = [];
    results.subject = [];
    results.cvholdout = [];
    results.finalholdout = [];
    results.lambda = [];
    results.lambda1 = [];
    results.LambdaSeq = [];
    results.Gtype = [];
    results.bias = [];
    results.normalize = [];
    results.data_varname = [];
    results.filters = [];
    results.p2 = [];
    results.cor2 = [];
    results.p2t = [];
    results.cor2t = [];
    results.err2 = [];
end
