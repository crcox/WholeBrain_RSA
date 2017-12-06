function [results,AdlasInstances] = learn_similarity_encoding(AdlasInstances, C, V, regularization, varargin)
% TO DO:
% Currently, the function can only handle a single subject. It will also
% only support HYPERBAND for a single hyperparameter, meaning that it can
% only do group lasso. This is because the function was originally written
% to facilitate grid search, and with HYPERBAND the hyperparameters should
% be be "crossed" in the same way. With HYPERBAND, configurations can be
% crossed with subject and cvholdout, and probably permutation index.

    p = inputParser();
    addRequired(p  , 'AdlasInstances');
    addRequired(p  , 'C');
    addRequired(p  , 'V');
    addRequired(p  , 'regularization');
    addParameter(p , 'cvind'                   , []       );
    addParameter(p , 'permutations'            , []       );
    addParameter(p , 'AdlasOpts'               , struct() );
    addParameter(p , 'Verbose'                 , true     );
    parse(p, C, V, regularization, target_type, varargin{:});

    if numel(p.Results.V) > 1
        error('crcox:TooManySubjects', 'Each call to "learn_similarity_encoding" currently supports only a single subject.')
    end
    if numel(p.Results.lambda) == numel(p.Results.lambda1) && numel(p.Results.lambda) > 1
        warning('crcox:NotImplemented', 'Note that HYPERBAND is not yet implemented for GrOWL. Lambda and Lambda1 will be crossed, as if grid searching.');
    end
    
    AdlasInstances = p.Results.AdlasInstances;
    C              = p.Results.C;
    V              = p.Results.V;
    regularization = p.Results.regularization;
    cvind          = p.Results.cvind;
    permutations   = p.Results.permutations;
    options        = p.Results.AdlasOpts;
    VERBOSE        = p.Results.Verbose;
    

    if strcmp(regularization, {'grOWL','grOWL2'});
        assert(~isempty(LambdaSeq),'A LambdaSeq type (linear or exponential) must be set when using grOWL*');
    end

    Vorig = V;

    if VERBOSE
        fprintf('%5s%6s%11s %11s  %11s %11s  \n', 'cv','lam','lam1','test err','train err','n vox')
    end

    for i = 1:numel(AdlasInstances)
        subix = AdlasInstances(i).subject;
        permix = AdlasInstances(i).RandomSeed;
        cvix = AdlasInstances(i).cvholdout;
        regularization = AdlasInstances(i).regularization;
        normalize = AdlasInstances(i).normalize;
        BIAS = AdlasInstances(i).bias;
        permutation_index = permutations{subix}(:,permix);

        test_set  = cvind{subix} == cvix; % CHECK THIS
        train_set = ~test_set;

        V = Vorig{subix};
        switch lower(normalize)
            case 'zscore_train'
                mm = mean(V(train_set,:),1);
                ss = std(V(train_set,:),0,1);
            case 'zscore'
                mm = mean(V,1);
                ss = std(V,0,1);
            case 'stdev_train'
                mm = zeros(1, size(V,2));
                ss = std(V(train_set,:),0,1);
            case 'stdev'
                mm = zeros(1, size(V,2));
                ss = std(V,0,1);
            case '2norm_train'
                mm = mean(V(train_set,:),1);
                ss = norm(V(train_set,:));
            case '2norm'
                mm = mean(V,1);
                ss = norm(V);
            otherwise
                error('Unrecognized normalizaion method! Exiting...')
        end
        z = ss > 0;
        V(:,z) = bsxfun(@minus,V(:,z), mm(z));
        V(:,z) = bsxfun(@rdivide,V(:,z), ss(z));
        if any(~z)
            warning('There are %d constant-valued voxels. These voxels are not normalized.', sum(z));
            if VERBOSE
                fprintf('Constant-valued voxel indexes:\n');
                disp(find(~z));
            end
        end

        if BIAS
            V = [V, ones(size(V,1),1)]; %#ok<AGROW> It's not actually growing, see line 86.
        end
        [~,d] = size(V);

        switch upper(regularization)
            case 'L1L2'
                lamseq = lam;
            case {'GROWL','GROWL2'} % There is no real distinction between GROWL and GROWL2 anymore, but for continuity I'll make GROWL2 map to this anyway.
                switch lower(LambdaSeq)
                    case 'linear'
                        lamseq = lam1*(d:-1:1)/d + lam;
                    case 'exponential'
                        lamseq = lam*sqrt(2*log((d*ones(1,d))./(1:d)));
                    case 'inf' % needs both lam and lam1
                        lamseq = [lam+lam1, repmat(lam,1,d-1)];
                end
            otherwise
                error('%s is not an implemented regularization. check spelling', regularization);
        end
        
        Cp = C{subix}(permutation_index,:);
        % cpcr = zeros(1,size(Cp,2));
        % for k = 1:size(Cp, 2);
        %     cpcr(k) = corr(Cp(:,k),C{subix}(:,k));
        % end
        % disp(cpcr);
        Ct = Cp(train_set,:);
        Ch = Cp(test_set,:);
        Vt = V(train_set,:);
                
        if isempty(AdlasInstances(iii).Adlas)
        % In case of new model:
        %   1. Initialize model
        %   2. Train model
            AdlasInstances(iii).Adlas = Adlas(Vt, Ct, lamseq, options);
            AdlasInstances(iii).Adlas = AdlasInstances(iii).Adlas.train(options);
        elseif AdlasInstances(iii).status == 2
        % In case of existing model:
        %   1. Check that status == 2, which means that the previos round
        %   of training stopped because it hit the iteration limit.
        %   2. If so, train for more iterations.
            AdlasInstances(iii).Adlas = AdlasInstances(iii).Adlas.train(options);
        else
        % In case of existing model and status == 1 or status == 3,
        % continue without doing anything. Status 1 means optimal
        % convergence with at least one nonzero weight assigned. Status 3
        % means the solution is zero-sparse.
        %
        % If a zero-sparse solution is obtained in a single iteration, this
        % can prevent some information from ever being logged in the Adlas
        % structure, which results in an error at run-time when trying to
        % pick up where things left off.
        %
        % By not trying to train these models any more (which is fine,
        % because they have converged on a solution already, anyway), this
        % error should be avoided.
        end

        Uz = AdlasInstances(iii).Adlas.X;
        info = AdlasInstances(iii).Adlas;

        % KEY
        % ===
        % Our models assume that S = V*W*V'
        % V  : The fMRI data.
        % W  : A matrix of weights.
        % S  : The true similarity matrix.
        % Sz : The predicted similarity matrix.
        % C  : The square root of S, truncated to the first r columns (low rank assumption)
        % Cz : The predicted C.
        % St : The approximated S, reconstructed from actual C.
        % Uz : The estimated voxel weights, with a r weights per voxel.

        nz_rows = any(Uz,2);
        ix = find(nz_rows);
        nv = size(Uz,1);
        Unz = nnz(nz_rows);
        uz = Uz(ix,:);

        % Prepare to evaluate solutions
%         Wz = uz*uz';
%         Sz = V(:,ix)*Wz*V(:,ix)';
        Cz = V(:,ix)*uz;
        Ctz = Cz(train_set,:);
        Chz = Cz(test_set,:);

        if any(test_set)
            err1 = nrsa_loss(Ch, Chz);
            %cor1 = nrsa_corr(Ch*Ch', Chz*Chz');
        else
            err1 = [];
            %cor1 = [];
        end
        err2 = nrsa_loss(Ct, Ctz);
        %cor2 = nrsa_corr(Ct*Ct', Ctz*Ctz');

        if VERBOSE
            fprintf('%3d | %6.2f | %6.2f | %10.2f | %10.2f | %10d | %10d\n', ...
                icv,lam,lam1,err1,err2,Unz,nv);
        end
    end
    fprintf('Exit status -- %s (%d iterations)\n', info.message, info.iter);
end % learn_similarity_encoding

% OLD PERMUTATION ALGORITHM
% -------------------------
%     PERMUTATION_INDEXES = cell(1, max(cvind));
%     for ic = unique(cvind)'
%         disp(sprintf('Permuting CV %d...', ic));
%         c = C(cvind==ic,:);
%         n = size(c,1);
%         if VERBOSE
%             fprintf('Permuting %d rows of C, independently by its %d columns.\n', n, r);
%             fprintf('First 10 rows of C, before shuffling.\n')
%             disp(c)
%         end
%         permix = randperm(n);
%         C(cvind==ic, :) = c(permix, :);
%         if VERBOSE
%             fprintf('First 10 rows of C, after shuffling.\n')
%             disp(c(permix,:))
%         end
%     end

% COMPARISON AGAINST FULL RANK AND LIMITED RANK SQUARE MATRIX
% ===========================================================
% St = C{subix}*C{subix}';
% if strcmpi(target_type,'similarity')
%     lt1  = logical(tril(true(nnz(test_set)),0));
%     s1   = S{subix}(test_set,test_set);
%     sz1  = Sz(test_set,test_set);
%     st1  = St(test_set,test_set);
%     s1   = s1(lt1);
%     sz1  = sz1(lt1);
%     st1  = st1(lt1);
% 
%     lt2 = logical(tril(true(nnz(train_set)),0));
%     s2  = S{subix}(train_set,train_set);
%     sz2 = Sz(train_set,train_set);
%     st2 = St(train_set,train_set);
%     s2  = s2(lt2);
%     sz2 = sz2(lt2);
%     st2 = st2(lt2);
% end

% COLLECT RESULTS
% ===============
% if ~SMALL
%     results(iii).Uz  = uz;
%     results(iii).Uix = uint32(ix(:)');
%     results(iii).Cz  = Cz;
%     results(iii).Sz  = Sz;
% end
% % Metadata
% results(iii).nzv            = uint32(Unz); % number of nonzero rows
% results(iii).nvox           = uint32(nv); % total number of voxels
% results(iii).subject        = subix;
% results(iii).cvholdout      = icv; % cross validation index
% results(iii).finalholdout   = []; % handled in parent function
% if strcmpi(regularization, 'L1L2_GLMNET')
%     results(iii).lambda1    = info.lambda;
% else
%     results(iii).lambda1    = lam1;
% end
% results(iii).lambda         = lam;
% results(iii).LambdaSeq      = LambdaSeq;
% results(iii).regularization = regularization;
% results(iii).bias           = BIAS;
% results(iii).normalize      = normalize;
% 
% if any(test_set)
%     results(iii).err1 = norm(Ch - Chz,'fro')/norm(Ch,'fro');
% end
% results(iii).err2 = norm(Ct - Ctz,'fro')/norm(Ct,'fro');
% 
% if isempty(lambda)
%     lambda_j = nan;
% else
%     lambda_j = lambda(j);
% end
% if isempty(lambda1)
%     lambda1_k = nan;
%     if strcmpi(regularization, 'L1L2_GLMNET')
%         lambda1_k = info.lambda;
%     end
% else
%     lambda1_k = lambda1(k);
% end
% 
% if strcmpi(regularization, 'L1L2_GLMNET')
%     results(iii).iter = info.npasses;
%     info.iter = info.npasses;
% else
%     results(iii).iter = info.iter;
% end