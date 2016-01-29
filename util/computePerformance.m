function results = computePerformance(X, C, CV, results, metadata, dataname)
% Meant for recomputing performance values, perhaps by applying weights to a
% different dataset.
  dvn = results.data_varname;
  for i = 1:numel(results.filters);
    lab = results.filters{i};
    z = strcmp({metadata.filter.label}, lab) & ...
        strcmp({metadata.filter.subset}, dvn) & ...
        [metadata.filter.dimension] == 2;
    if any(z)
      zc = metadata.filter(z).filter;
      break
    end
  end

  for i = 1:numel(results.filters);
    lab = results.filters{i};
    z = strcmp({metadata.filter.label}, lab) & ...
        strcmp({metadata.filter.subset}, dataname) & ...
        [metadata.filter.dimension] == 1;
    if any(z)
      zr = metadata.filter(z).filter;
      break
    end
  end
  X = X(zr, zc);
  if results.bias
    X = [X,ones(size(X,1),1)];
  end

  Cz = X * results.Uz;
  C = C(zr,:);

  if all(size(CV)>1)
    % Passed matrix of CV schemes
    error('WholeBrain_RSA:computePerformance:invalidInput', ...
          '*** Error: CV must be a vector. ***');
  end

  CV = CV(zr);
  omit = CV == results.finalholdout;
  test = CV == results.cvholdout;
  train = ~test & ~omit;

  results.err1 = norm( C(test,:)-Cz(test,:) ,'fro') / norm( C(test,:) , 'fro');
  results.err2 = norm( C(train,:)-Cz(train,:) ,'fro') / norm( C(train,:) , 'fro');
end
