function [results,AdlasInstances] = learn_similarity_encoding(C, V, regularization, target_type, varargin)
% TO DO:
% Currently, the function can only handle a single subject. It will also
% only support HYPERBAND for a single hyperparameter, meaning that it can
% only do group lasso. This is because the function was originally written
% to facilitate grid search, and with HYPERBAND the hyperparameters should
% be be "crossed" in the same way. With HYPERBAND, configurations can be
% crossed with subject and cvholdout, and probably permutation index.

    % Local Constants:
    % ---------------
    GLMNET_ALPHA_LASSO = 1;
%     GLMNET_ALPHA_RIDGE = 0;

    p = inputParser();
    addRequired(p  , 'C');
    addRequired(p  , 'V');
    addRequired(p  , 'regularization');
    addRequired(p  , 'target_type');
    addParameter(p , 'tau'                     , []       );
    addParameter(p , 'lambda'                  , []       );
    addParameter(p , 'lambda1'                 , []       );
    addParameter(p , 'LambdaSeq'               , []       );
    addParameter(p , 'cvind'                   , []       );
    addParameter(p , 'cvholdout'               , []       );
    addParameter(p , 'normalize'               , []       );
    addParameter(p , 'bias'                    , []       );
    addParameter(p , 'permutations'            , []       );
    addParameter(p , 'DEBUG'                   , false    );
    addParameter(p , 'AdlasOpts'               , struct() );
    addParameter(p , 'SmallFootprint'          , false    );
    addParameter(p , 'Verbose'                 , true     );
    addParameter(p , 'PARALLEL'                , false    );
    addParameter(p , 'AdlasInstances'          , []       );
    parse(p, C, V, regularization, target_type, varargin{:});

    if numel(p.Results.V) > 1
        error('crcox:TooManySubjects', 'Each call to "learn_similarity_encoding" currently supports only a single subject.')
    end
    if numel(p.Results.lambda) == numel(p.Results.lambda1) && numel(p.Results.lambda) > 1
        warning('crcox:NotImplemented', 'Note that HYPERBAND is not yet implemented for GrOWL. Lambda and Lambda1 will be crossed, as if grid searching.');
    end
    C                       = p.Results.C;
    V                       = p.Results.V;
    regularization          = p.Results.regularization;
    lambda                  = p.Results.lambda;
    lambda1                 = p.Results.lambda1;
    LambdaSeq               = p.Results.LambdaSeq;
    cvind                   = p.Results.cvind;
    holdout                 = p.Results.cvholdout;
    normalize               = p.Results.normalize;
    permutations            = p.Results.permutations;
    BIAS                    = p.Results.bias;
    DEBUG                   = p.Results.DEBUG;
    options                 = p.Results.AdlasOpts;
    SMALL                   = p.Results.SmallFootprint;
    VERBOSE                 = p.Results.Verbose;
    PARALLEL                = p.Results.PARALLEL;

    if strcmp(regularization, {'grOWL','grOWL2'});
        assert(~isempty(LambdaSeq),'A LambdaSeq type (linear or exponential) must be set when using grOWL*');
    end

    Vorig = V;

    if isempty(lambda)
        nlam = 1;
    else
        nlam = length(lambda);
    end

    if isempty(lambda1)
        nlam1 = 1;
    else
        nlam1 = length(lambda1);
    end

    if isempty(holdout)
        cvset = 1:max(cvind);
    else
        cvset = holdout;
    end
    ncv = numel(cvset);

    % Define results structure
    results = struct( ...
        'Uz'             , [] , ...
        'Cz'             , [] , ...
        'Sz'             , [] , ...
        'nz_rows'        , [] , ...
        'subject'        , [] , ...
        'cvholdout'      , [] , ...
        'finalholdout'   , [] , ...
        'lambda'         , [] , ...
        'lambda1'        , [] , ...
        'LambdaSeq'      , [] , ...
        'regularization' , [] , ...
        'bias'           , [] , ...
        'normalize'      , [] , ...
        'nzv'            , [] , ...
        'coords'         , [] , ...
        'err1'           , [] , ...
        'err2'           , [] , ...
        'iter'           , [] );

    % Preallocate
    if isempty(permutations)
        nperm = 1;
    else
        nperm = size(permutations{1}, 2);
    end
    N = numel(cvset)*nlam*nperm;
    results(N).Uz = [];

    if isempty(p.Results.AdlasInstances)
        AdlasInstances = repmat(Adlas,N,1);
    else
        AdlasInstances = p.Results.AdlasInstances;
    end

    if VERBOSE
        fprintf('%5s%6s%11s %11s  %11s %11s  \n', 'cv','lam','lam1','test err','train err','n vox')
    end

    iii = 0; % index into 1-D results structure.
    for subix = 1:numel(Vorig)
        for permix = 1:nperm
            permutation_index = permutations{subix}(:,permix);
            for i = 1:ncv
                icv = cvset(i);
                test_set  = cvind{subix}==icv;
                train_set = ~test_set;

                V = Vorig{subix};

                % normalize
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
                    V = [V, ones(size(V,1),1)]; %#ok<AGROW> It's not actually growing, see line 111.
                end

                [~,d] = size(V);

                Cp = C{subix}(permutation_index,:);
                % cpcr = zeros(1,size(Cp,2));
                % for k = 1:size(Cp, 2);
                %     cpcr(k) = corr(Cp(:,k),C{subix}(:,k));
                % end
                % disp(cpcr);
                Ct = Cp(train_set,:);
                Ch = Cp(test_set,:);
                Vt = V(train_set,:);
                Vh = V(test_set,:);

                for j = 1:nlam
                    if isempty(lambda)
                        lam = nan(1);
                    else
                        lam = lambda(j);
                    end

                    for k = 1:nlam1
                        iii = iii + 1;
                        if isempty(lambda1)
                            lam1 = nan(1);
                        else
                            lam1 = lambda1(k);
                        end

                        if DEBUG
                            Uz = randn(d,r);
                            info.message = 'DEBUG';
                            info.iter    = 0;
                        else
                            switch upper(regularization)
                                case 'L1L2_GLMNET'
                                    cvind_glmnet = drop_index_and_adjust(cvind{subix}, icv);
                                    [Uz, info] = grouplasso_glmnet(Vt, Ct, ...
                                        GLMNET_ALPHA_LASSO, ...
                                        lam,...
                                        'cvind'   , cvind_glmnet,...
                                        'bias'    , BIAS, ...
                                        'PARALLEL', PARALLEL, ...
                                        'tol'     , 1e-8, ...
                                        'U0'      ,   []);

                                case 'L1L2'
                                    if all(lam==0)
                                        Uz = pinv(Vt)*Ct;
                                        info = struct();
                                    else
                                        if isempty(AdlasInstances(iii))
                                            AdlasInstances(iii) = Adlas(Vt, Ct, lam, options);
                                        end
                                        AdlasInstances(iii) = AdlasInstances(iii).train(options);
                                        Uz = AdlasInstances(iii).X;
                                        info = AdlasInstances(iii);
                                    end

                                case 'GROWL'
                                    switch LambdaSeq
                                        case 'linear'
                                            lamseq = lam*(d:-1:1)/d + lam1;
                                        case 'exponential'
                                            lamseq = lam*sqrt(2*log((d*ones(1,d))./(1:d)));
                                    end
                                    if all(lamseq==0)
                                        Uz = pinv(Vt)*Ct;
                                        info = struct();
                                    else
                                        [Uz, info] = Adlas1(Vt, Ct, lamseq, options);
                                    end

                                case 'GROWL2'
                                    switch LambdaSeq
                                        case 'linear'
                                            lamseq = lam*(d:-1:1)/d;
                                        case 'exponential'
                                            lamseq = lam*sqrt(2*log((d*ones(1,d))./(1:d)));
                                        case 'inf' % needs both lam and lam1
                                            lamseq = [lam+lam1, repmat(lam,1,d-1)];
                                    end
                                    if all(lamseq==0)
                                        Uz = pinv(Vt)*Ct;
                                        info = struct();
                                    else
                                        [Uz, info] = Adlas1(Vt, Ct, lamseq, options);
                                    end
                                otherwise
                                    error('%s is not an implemented regularization. check spelling', regularization);
                            end
                        end

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
                        Wz = uz*uz';
                        Sz = V(:,ix)*Wz*V(:,ix)';
                        Cz = V(:,ix)*uz;
                        Ctz = Vt(:,ix)*uz;
                        Chz = Vh(:,ix)*uz;
    %                     St = C{subix}*C{subix}';

    %                     if strcmpi(target_type,'similarity')
    %                         lt1  = logical(tril(true(nnz(test_set)),0));
    %                         s1   = S{subix}(test_set,test_set);
    %                         sz1  = Sz(test_set,test_set);
    %                         st1  = St(test_set,test_set);
    %                         s1   = s1(lt1);
    %                         sz1  = sz1(lt1);
    %                         st1  = st1(lt1);
    % 
    %                         lt2 = logical(tril(true(nnz(train_set)),0));
    %                         s2  = S{subix}(train_set,train_set);
    %                         sz2 = Sz(train_set,train_set);
    %                         st2 = St(train_set,train_set);
    %                         s2  = s2(lt2);
    %                         sz2 = sz2(lt2);
    %                         st2 = st2(lt2);
    %                     end

                        if ~SMALL
                            results(iii).Uz  = uz;
                            results(iii).Uix = uint32(ix(:)');
                            results(iii).Cz  = Cz;
                            results(iii).Sz  = Sz;
                        end
                        % Metadata
                        results(iii).nzv            = uint32(Unz); % number of nonzero rows
                        results(iii).nvox           = uint32(nv); % total number of voxels
                        results(iii).subject        = subix;
                        results(iii).cvholdout      = icv; % cross validation index
                        results(iii).finalholdout   = []; % handled in parent function
                        if strcmpi(regularization, 'L1L2_GLMNET')
                            results(iii).lambda1    = info.lambda;
                        else
                            results(iii).lambda1    = lam1;
                        end
                        results(iii).lambda         = lam;
                        results(iii).LambdaSeq      = LambdaSeq;
                        results(iii).regularization = regularization;
                        results(iii).bias           = BIAS;
                        results(iii).normalize      = normalize;

                        if any(test_set)
                            results(iii).err1 = norm(Ch - Chz,'fro')/norm(Ch,'fro');
    %                         if strcmpi(target_type,'similarity')
    %                             results(iii).structureScoreMap = s1(:)' * sz1(:);
    %                         end
                        end

                        results(iii).err2 = norm(Ct - Ctz,'fro')/norm(Ct,'fro');
    %                     if strcmpi(target_type,'similarity')
    %                         results(iii).structureScoreMap = s2(:)' * sz2(:);
    %                     end

                        if strcmpi(regularization, 'L1L2_GLMNET')
                            results(iii).iter = info.npasses;
                            info.iter = info.npasses;
                        else
                            results(iii).iter = info.iter;
                        end

                        if isempty(lambda)
                            lambda_j = nan;
                        else
                            lambda_j = lambda(j);
                        end
                        if isempty(lambda1)
                            lambda1_k = nan;
                            if strcmpi(regularization, 'L1L2_GLMNET')
                                lambda1_k = info.lambda;
                            end
                        else
                            lambda1_k = lambda1(k);
                        end

                        err1 = results(iii).err1;
                        err2 = results(iii).err2;
                        if VERBOSE
                            fprintf('%3d | %6.2f | %6.2f | %10.2f | %10.2f | %10d\n', ...
                                icv, lambda_j,lambda1_k,err1,err2,Unz);
                        end
                    end % lam1 loop
                end % lam loop
            end % cv loop
        end
    end % subject loop
    if VERBOSE
        fprintf('Exit status -- %s (%d iterations)\n', info.message, info.iter);
    end
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
