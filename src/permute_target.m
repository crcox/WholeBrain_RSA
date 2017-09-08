function X = permute_target(C,method,arg)
% PERMUTE_TARGET Needs documentation
%
% Chris Cox 25/08/2017
    if ~iscell(C)
        C = {C};
    end
    if nargin < 2
        method = 'simple';
        cvind = cell(size(C));
        for i = 1:numel(C);
            cvind{i} = ones(size(C{i},1),1);
        end
        permix = [];
    else
        switch method
            case {'simple', 'simulated'}
                if nargin > 2
                    error('crcox:ToManyArgs', 'The %s method does not accept any arguments.', method);
                end
                cvind = cell(size(C));
                for i = 1:numel(C);
                    cvind{i} = ones(size(C{i},1),1);
                end
                permix = [];
            case 'stratified'
                if nargin < 3
                    error('crcox:ToFewArgs', 'The %s method requires an argument that specifies the cross validation structure.', method);
                end
                if iscell(arg);
                    cvind = arg;
                else
                    cvind = {arg};
                end
                permix = [];
            case 'manual'
                if nargin < 3
                    error('crcox:ToFewArgs', 'The %s method requires an argument that specifies the permutation index.', method);
                end
                cvind = cell(size(C));
                for i = 1:numel(C);
                    cvind{i} = ones(size(C{i},1),1);
                end
                permix = arg;
        end
    end 

    X = cell(size(C));
    for i = 1:numel(C)
        if size(C{i},1) == 1
            C{i} = C{i}';
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
                case 'manual'
                    X{i} = manual(y, permix);
                    break % Don't need the CV loop for this method
            end
        end
    end
    if iscell(X) && numel(X) == 1 && ~iscell(C);
        X = X{1};
    end
end

function x = simple_perm(y)
% Takes y and randomly permutes it
    permix = randperm(size(y,1));
    x = y(permix, :);
end

function [x,ix] = stratified_perm(y, g)
% Stratified permuation involves permuting within each category, so that
% each category member is equally likely to be reassigned to each other
% group. Another way to put it is that the permuted result will be as close
% to orthogonal with the true structure as possible, which means that the
% hypothesized structure between the data and this target structure will
% have been eliminated.

% This is not implemented (yet?) because it make less intuitive sense when
% the targets are continuous and not fundamentally categorical.
    warning('stratified permutation is not implemented yet. Doing simple permutation instead...');
    [x,ix] = simple_perm(y);
end

function [x,ix] = simulated(y)
% Simulates a new y based on its mean and covariance.
    x = mvnrnd(mean(y), cov(y), size(y,1));
    ix = nan(1);
end

function [x,ix] = manual(y, permix)
    ix = permix;
    x = y(ix, :);
end