function Xa = averageRepeatedTrials(X, stimcode)
%AVERAGEREPEATEDTRIALS  Define k-fold or LOO cross validation assignments
%
% averageRepeatedTrials(X, stimcode)
%
% stimcode is either a numeric vector or a cell array of strings that
% identify each stimulus. This function will average rows that share the
% same stimcode.
% 
% Examples:
%   X = randn(16,10);
%   stimcode = repmat(1:8,1,2);
%   averageRepeatedTrials(X, stimcode)
%
% See also:
% BAR
% SOMECLASS/SOMEMETHOD

  % Check that X and stimcode are compatible
  if size(X,1) ~= numel(stimcode)
    err.message = sprintf('size(X,1) must equal numel(stimcode).');
    err.identifier = 'WholeBrain_RSA:averageRepeatedTrials:incompatibleArguments';
    error(err);
  end
  
  % Parse stimcode and generate stimid (internal, integer representation).
  stimset = unique(stimcode, 'sorted');
  stimid = zeros(numel(stimcode),1);
  for i = 1:numel(stimset)
    if iscellstr(stimset)
      z = strcmp(stimset{i}, stimcode);
    else %isnumeric
      z = stimset(i) == stimcode;
    end
    stimid(z) = i;
  end
  
  % Perform averaging
  n = max(stimid);
  Xa = zeros(n, size(X,2));
  for i = 1:n
    z = (stimid == i);
    Xa(i,:) = mean(X(z,:));
  end
end
    
  