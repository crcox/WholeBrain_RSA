function [sigvox] = ThresholdNodestrength(nodestrength, empiricalnullvals, varargin)
  p = inputParser();
  addRequired(p, 'nodestrength', @isnumeric)
  addRequired(p, 'empiricalnullvals', @isnumeric)
  addParameter(p, 'method', 'bonferroni', @ischar)
  addParameter(p, 'alpha', 0.05, @isnumeric)
  parse(p, nodestrength, empiricalnullvals, varargin{:});

  nodestrength = p.Results.nodestrength;
  empiricalnullvals = p.Results.empiricalnullvals;
  method = lower(p.Results.method);
  alpha = lower(p.Results.alpha);

  n = numel(nodestrength);
  %n = nnz(nodestrength);
  switch method
  case 'fdr'
    % Threshold the voxel with the largest node strength by (1*alpha)/n, the
    % voxel with the second largest node strength with (2*alpha)/n, and so on
    % until a voxel fails to pass threshold.
    [nodestrength_sorted, ix] = sort(nodestrength,'descend');
    [~,ix_unsort] = sort(ix);
    empiricalnullvals_sorted = empiricalnullvals(:,ix);
    alpha_fdr = ((1:n) / n) * alpha;
    sigvox = false(1,n);
    for i = 1:n;
      thresh = quantile(empiricalnullvals_sorted(:,i), 1-alpha_fdr(1));
      fprintf('%d, %.8f, %.8f, %.8f\n', i, 1-alpha_fdr(i), thresh, nodestrength_sorted(i));
      if nodestrength_sorted(i) > thresh
        sigvox(i) = true;
      else
        break;
      end
    end
    sigvox = sigvox(ix_unsort);

  case 'bonferroni'
    % Threshold all voxels by alpha/n.
    alpha_bonf = alpha / n;
    disp(alpha_bonf);
    thresh = quantile(empiricalnullvals, 1-alpha_bonf);
    sigvox = nodestrength(:) > thresh(:);

  otherwise
    printf('*** Unknown method "%s". Choose fdr or bonferroni. ***\n');
    error('WholeBrain_RSA:ThresholdNodestrength:unknownmethod');

end
