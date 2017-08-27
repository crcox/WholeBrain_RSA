function [results,info] = learn_similarity_encoding(S, V, regularization, target_type, varargin)
    % Local Constants:
    % ---------------
    GLMNET_ALPHA_LASSO = 1;
    GLMNET_ALPHA_RIDGE = 0;

    p = inputParser();
    addRequired(p  , 'S');
    addRequired(p  , 'V');
    addRequired(p  , 'regularization');
    addRequired(p  , 'target_type');
    addParameter(p , 'tau'                     , []       );
    addParameter(p , 'lambda'                  , []       );
    addParameter(p , 'lambda1'                 , []       );
    addParameter(p , 'cvind'                   , []       );
    addParameter(p , 'cvholdout'               , []       );
    addParameter(p , 'normalize'               , []       );
    addParameter(p , 'bias'                    , []       );
    addParameter(p , 'LambdaSeq'               , []       );
    addParameter(p , 'DEBUG'                   , false    );
    addParameter(p , 'PermutationTest'         , false    );
    addParameter(p , 'PermutationMethod'       , 'simple' );
    addParameter(p , 'RestrictPermutationByCV' , false    );
    addParameter(p , 'AdlasOpts'               , struct() );
    addParameter(p , 'SmallFootprint'          , false    );
    addParameter(p , 'Verbose'                 , true     );
    addParameter(p , 'PARALLEL'                , false    );
    parse(p, S, V, regularization, target_type, varargin{:});

    S                       = p.Results.S;
    V                       = p.Results.V;
    regularization          = p.Results.regularization;
    target_type             = p.Results.target_type;
    tau                     = p.Results.tau;
    lambda                  = p.Results.lambda;
    lambda1                 = p.Results.lambda1;
    cvind                   = p.Results.cvind;
    holdout                 = p.Results.cvholdout;
    normalize               = p.Results.normalize;
    BIAS                    = p.Results.bias;
    LambdaSeq               = p.Results.LambdaSeq;
    PermutationTest         = p.Results.PermutationTest;
    PermutationMethod       = p.Results.PermutationMethod;
    RestrictPermutationByCV = p.Results.RestrictPermutationByCV;
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
    results(numel(cvset)*nlam1*nlam).Uz = [];

    %square root
    for i = 1:numel(S)
        switch target_type
            case 'similarity'
                [C{i}, r] = sqrt_truncate_r(S{i}, tau);
                if VERBOSE
                    fprintf('S decomposed into %d dimensions (tau=%.2f)\n', r, tau)
                end
            case 'embedding'
                C{i} = S{i};
                r = size(C{i},2);
        end
    end

    if VERBOSE
        fprintf('PermutationTest: %d\n', PermutationTest);
    end
    if PermutationTest
        % Note on randomization for permutation
        % -------------------------------------
        % A required argument when specifying permutations is a list of
        % "RandomSeeds". These are applied near the beginning of the
        % program (within WholeBrain_RSA), to seed the random number
        % generator.
        %
        % If the PermutationMethod is 'manual', then the RandomSeed has a
        % different (additional) function. It will be used to index into
        % the columns of a n x p matrix, generated in advance, that
        % contains the indexes to generate p unique permutations.
        %
        % In this case, the matrix should stored in a variable named
        % PERMUTATION_INDEXES, contained within a file named
        % PERMUTATION_INDEXES.mat
        fprintf('PermutationTest: %d\n', PermutationTest);
        if PermutationTest
            switch PermutationMethod
                case 'simple'
                    if RestrictPermutationByCV
                        C = permute_target(C, PermutationMethod, cvind);
                    else
                        C = permute_target(C, PermutationMethod);
                    end
                case 'manual'
                    load('PERMUTATION_INDEXES.mat', 'PERMUTATION_INDEXES');
                    pix = PERMUTATION_INDEXES(:, RandomSeed);
                    C = permute_target(C, PermutationMethod, pix);
            end
        end
    end

    if VERBOSE
        fprintf('%5s%6s%11s %11s  %11s %11s  \n', 'cv','lam','lam1','test err','train err','n vox')
    end

    iii = 0; % index into 1-D results structure.
    for subix = 1:numel(Vorig)
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
                V = [V, ones(size(V,1),1)];
            end

            [~,d] = size(V);

            Ct = C{subix}(train_set,:);
            Ch = C{subix}(test_set,:);
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
                                    [Uz, info] = Adlas1(Vt, Ct, lam, options);
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
                    St = C{subix}*C{subix}';

                    if strcmpi(target_type,'similarity')
                        lt1  = logical(tril(true(nnz(test_set)),0));
                        s1   = S{subix}(test_set,test_set);
                        sz1  = Sz(test_set,test_set);
                        st1  = St(test_set,test_set);
                        s1   = s1(lt1);
                        sz1  = sz1(lt1);
                        st1  = st1(lt1);

                        lt2 = logical(tril(true(nnz(train_set)),0));
                        s2  = S{subix}(train_set,train_set);
                        sz2 = Sz(train_set,train_set);
                        st2 = St(train_set,train_set);
                        s2  = s2(lt2);
                        sz2 = sz2(lt2);
                        st2 = st2(lt2);
                    end

                    if ~SMALL
                        results(iii).Uz  = uz;
                        results(iii).Uix = uint32(ix(:)');
                        results(iii).Cz  = Cz;
                        results(iii).Sz  = Sz;
                    end
                    % Metadata
                    results(iii).nzv            = uint32(Unz); % number of nonzero rows
                    results(iii).nvox           = uint32(nv); % total number of voxels
                    results(iii).subject        = []; % handled in parent function
                    results(iii).cvholdout      = icv; % cross validation index
                    results(iii).finalholdout   = []; % handled in parent function
                    if strcmp(upper(regularization), 'L1L2_GLMNET')
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
                        if strcmpi(target_type,'similarity')
                            results(iii).structureScoreMap = s1(:)' * sz1(:);
                        end
                    end

                    results(iii).err2 = norm(Ct - Ctz,'fro')/norm(Ct,'fro');
                    if strcmpi(target_type,'similarity')
                        results(iii).structureScoreMap = s2(:)' * sz2(:);
                    end

                    if strcmp(upper(regularization), 'L1L2_GLMNET')
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
                        if strcmp(upper(regularization), 'L1L2_GLMNET')
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
