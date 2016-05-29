function Wz = computeWz(Uz,varargin)
  p = inputParser();
  addRequired(p, 'Uz', @ismatrix);
  addParameter(p, 'full', false, @islogical);
  parse(p,varargin{:});

  Uz = p.Results.Uz;
  IncludeZeroValues = p.Results.full;

  if IncludeZeroValues
    Wz = Uz * Uz';
  else
    z = any(Uz,2);
    Wz = Uz(z,:) * Uz(z,:)';
  end
end
