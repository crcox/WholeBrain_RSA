function X = permute_target(C,method,cvind)
    if ~iscell(C)
        C = {C};
    end
    if nargin < 2
        method = 'simple';
    end 
    if nargin < 3
        cvind = cell(size(C));
    else
        if ~iscell(cvind)
            cvind = {cvind};
        end
    end
    X = cell(size(C));
    for i = 1:numel(C)
        if size(C{i},1) == 1
            C{i} = C{i}';
        end
        if isempty(cvind{i})
            cvind{i} = ones(size(C{i},1),1);
        end
%         if strcmp('stratified', method)
%             % Done here so it's not recomputed for each holdout.
%             g = cluster(linkage(pdist(C{i},'cosine')), 'maxclust', 2);
%         end
        cvset = unique(cvind{i})';
        X{i} = zeros(size(C{i}));
        for j = 1:numel(cvset); % unique(cvind{i})'
            ic = cvset(j);
            y = C{i}(cvind{i}==ic,:);
            switch method
                case 'simple'
                    X{i}(cvind{i}==ic,:) = simple_perm(y);
                case 'stratified'
                    X{i}(cvind{i}==ic,:) = stratified_perm(y, g);
                case 'simulated'
                    X{i}(cvind{i}==ic,:) = simulated(y);
            end
        end
    end
    if iscell(X) && numel(X) == 1
        X = X{1};
    end
end

function x = simple_perm(y)
% Takes y and randomly permutes it
    permix = randperm(size(y,1));
    x = y(permix, :);
end

function x = stratified_perm(y, g)
% Stratified permuation involves permuting within each category, so that
% each category member is equally likely to be reassigned to each other
% group. Another way to put it is that the permuted result will be as close
% to orthogonal with the true structure as possible, which means that the
% hypothesized structure between the data and this target structure will
% have been eliminated.

% This is not implemented (yet?) because it make less intuitive sense when
% the targets are continuous and not fundamentally categorical.
    warning('stratified permutation is not implemented yet. Doing simple permutation instead...');
    x = simple_perm(y);
end

function x = simulated(y)
% Simulates a new y based on its mean and covariance.
    x = mvnrnd(mean(y), cov(y), size(y,1));
end