function response = QueryResults(results,params,varargin)
  p = inputParser();
  addParameter(p,'ModelType',@isValidModelType);
  addParameter(p,'DataType',@isValidDataType);
  addParameter(p,'Biased',@islogical);
  addParameter(p,'Normalized',@islogical);
  parse(p,varargin{:});

  mt = p.Results.ModelType;
  dt = p.Results.DataType;
  bias = p.Results.Biased;
  norm = p.Results.Normalized;

  z = strcmp(mt,{params.SanityCheckModel}) & ...
      strcmp(dt,{params.SanityCheckData}) & ...
      [params.normalize] == norm & ...
      [params.bias] == bias;

  response = results(z);
end

function isvalid = isValidModelType(type)
  ModelTypes = {'use_random_cor','use_random_inner','use_shuffled_cor','use_shuffled_inner','real'};
  isvalid = any(strcmp(type,ModelTypes));
  if ~isvalid
      error('%s is not a valid model type.\n',type)
  end
end
function isvalid = isValidDataType(type)
  DataTypes = {'userandom','useshuffled','real'};
  isvalid = any(strcmp(type,DataTypes));
  if ~isvalid
      error('%s is not a valid data type.\n',type)
  end
end
