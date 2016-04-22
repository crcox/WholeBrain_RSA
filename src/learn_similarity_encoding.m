function [results,info] = learn_similarity_encoding(S, V, regularization, varargin)
  p = inputParser();
  addRequired(p,'S');
  addRequired(p,'V');
  addRequired(p,'regularization');
  addParameter(p , 'tau'       , 0.2);
  addParameter(p , 'lambda'    , []);
  addParameter(p , 'lambda1'   , []);
  addParameter(p , 'cvind'     , []);
  addParameter(p , 'cvholdout' , []);
  addParameter(p , 'normalize' , []);
  addParameter(p , 'LambdaSeq' , []);
  addParameter(p , 'DEBUG'     , false);
  addParameter(p , 'PermutationTest'     , false);
  addParameter(p , 'AdlasOpts' , struct());
  addParameter(p , 'SmallFootprint' ,false);
  parse(p, S, V, regularization, varargin{:});

  S         = p.Results.S;
  V         = p.Results.V;
  regularization     = p.Results.regularization;
  tau       = p.Results.tau;
  lambda    = p.Results.lambda;
  lambda1   = p.Results.lambda1;
  cvind     = p.Results.cvind;
  holdout   = p.Results.cvholdout;
  normalize = p.Results.normalize;
  LambdaSeq = p.Results.LambdaSeq;
  PermutationTest = p.Results.PermutationTest;
  DEBUG     = p.Results.DEBUG;
  options   = p.Results.AdlasOpts;
  SMALL     = p.Results.SmallFootprint;

  if strcmp(regularization, {'grOWL','grOWL2'});
    assert(~isempty(LambdaSeq),'A LambdaSeq type (linear or exponential) must be set when using grOWL*');
  end

  [~,d] = size(V);
  Vorig = V;
  BIASUNIT = all(V(:,end)==1);

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
  results.Uz = [];
  results.Cz = [];
  results.Sz = [];
  results.nz_rows =  [];
  results.subject =  [];
  results.cvholdout = [];
  results.finalholdout = [];
  results.lambda = [];
  results.lambda1 = [];
  results.LambdaSeq = [];
  results.regularization = [];
  results.bias = [];
  results.normalize = [];
  results.nzv = [];
  results.p1      =  [];
  results.p2      =  [];
  results.cor1    =  [];
  results.cor2    =  [];
  results.p1t     =  [];
  results.p2t     =  [];
  results.cor1t   =  [];
  results.cor2t   =  [];
  results.err1    =  [];
  results.err2    =  [];
  results.iter    =  [];

  % Preallocate
  results(numel(cvset)*nlam1*nlam).Uz = [];

  %square root
  [C, r] = sqrt_truncate_r(S, tau);

  fprintf('PermutationTest: %d\n', PermutationTest);
  if PermutationTest
    n = size(C,1);
    fprintf('Permuting %d rows of C, independently by its %d columns.\n', n, r);
    fprintf('First 10 rows of C, before shuffling.\n')
    disp(C(1:10,:))
    permix = randperm(n);
    C = C(permix, :);
%    for i = 1:r
%      permix = randperm(n);
%      C(:,i) = C(permix, i);
%    end
    fprintf('First 10 rows of C, after shuffling.\n')
    disp(C(1:10,:))
  end

  fprintf('%5s%6s%11s %11s  %11s %11s  \n', 'cv','lam','lam1','test err','train err','n vox')

  iii = 0; % index into 1-D results structure.
  for i = 1:ncv
    icv = cvset(i);
    test_set  = cvind==icv;
    train_set = ~test_set;

    V = Vorig;

    % normalize
    switch lower(normalize)
    case 'zscore_train'
      mm = mean(V(train_set,1:nv),1);
      ss = std(V(train_set,1:nv),0,1);
    case 'zscore'
      mm = mean(V(:,1),1);
      ss = std(V,0,1);
    case 'stdev'
      mm = 0;
      ss = std(V,0,1);
    case '2norm'
      mm = mean(V,1);
      ss = norm(V);
    otherwise
      mm = 0;
      ss = 1;
    end
    V = bsxfun(@minus,V, mm);
    z = ss > 0;
    V(:,z) = bsxfun(@rdivide,V(:,z), ss(z));

    C1 = C(train_set,:);
    V1 = V(train_set,:);

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
          case 'L1L2'
            if all(lam==0)
              Uz = pinv(V1)*C1;
              info = struct();
            else
              [Uz, info] = Adlas1(V1, C1, lam1, options);
            end

          case 'GROWL'
            switch LambdaSeq
            case 'linear'
              lamseq = lam*(d:-1:1)/d;
            case 'exponential'
              lamseq = lam*sqrt(2*log((d*ones(1,d))./(1:d)));
            end
            if all(lamseq==0)
              Uz = pinv(V1)*C1;
              info = struct();
            else
              [Uz, info] = Adlas1(V1, C1, lamseq, options);
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
              Uz = pinv(V1)*C1;
              info = struct();
            else
              [Uz, info] = Adlas1(V1, C1, lamseq, options);
            end
          otherwise
            error('%s is not an implemented regularization. check spelling', regularization);
          end
        end

        % KEY:
        % Our models assume that S = V*W*V'
        % V  : The fMRI data.
        % W  : A matrix of weights.
        % S  : The true similarity matrix.
        % Sz : The predicted similarity matrix.
        % C  : The square root of S, truncated to the first r columns (low rank assumption)
        % Cz : The predicted C.
        % St : The approximated S, reconstructed from actual C.
        % Uz : The estimated voxel weights, with a r weights per voxel.

        % Prepare to evaluate solutions
        k1 = nnz(any(Uz,2));
        Wz = Uz*Uz';
        Sz = V*Wz*V';
        Cz = V*Uz;
        St = C*C';

        lt1  = logical(tril(true(nnz(test_set)),0));
        s1   = S(test_set,test_set);
        sz1  = Sz(test_set,test_set);
        st1  = St(test_set,test_set);
        s1   = s1(lt1);
        sz1  = sz1(lt1);
        st1  = st1(lt1);

        lt2 = logical(tril(true(nnz(train_set)),0));
        s2  = S(train_set,train_set);
        sz2 = Sz(train_set,train_set);
        st2 = St(train_set,train_set);
        s2  = s2(lt2);
        sz2 = sz2(lt2);
        st2 = st2(lt2);

        if ~SMALL
          results(iii).Uz = Uz;
          results(iii).Cz = Cz;
          results(iii).Sz = Sz;
          results(iii).nz_rows = any(Uz,2);
        end
        % Metadata
        results(iii).subject =  []; % handled in parent function
        results(iii).cvholdout = icv;
        results(iii).finalholdout = []; % handled in parent function
        results(iii).lambda = lam;
        results(iii).lambda1 = lam1;
        results(iii).LambdaSeq = LambdaSeq;
        results(iii).regularization = regularization;
        results(iii).bias = BIASUNIT;
        results(iii).normalize = normalize;
        results(iii).nzv = k1;

        if any(test_set)
%          results(iii).p1      = trace(corr(S(test_set,:)',Sz(test_set,:)'))/nnz(test_set);
%          results(iii).cor1    = corr(s,sz); % test
%          results(iii).p1t     = trace(corr(St(test_set,:)',Sz(test_set,:)'))/nnz(test_set);
%          results(iii).cor1t   = corr(st,sz); % test
          results(iii).err1    = norm(C(test_set,:) - Cz(test_set,:),'fro')/norm(C(test_set,:),'fro');
          results(iii).structureScoreMap = s1(:)' * sz1(:);
        end

%        results(iii).p2      = trace(corr(S(train_set,:)',Sz(train_set,:)'))/nnz(train_set);
%        results(iii).cor2    = corr(s2,sz2); % train
%        results(iii).p2t     = trace(corr(St(train_set,:)',Sz(train_set,:)'))/nnz(train_set);
%        results(iii).cor2t   = corr(st2,sz2); % train
        results(iii).err2    = norm(C(train_set,:) - Cz(train_set,:),'fro')/norm(C(train_set,:),'fro');
        results(iii).structureScoreMap = s2(:)' * sz2(:);

        results(iii).iter    = info.iter;

        if isempty(lambda)
          lambda_j = nan;
        else
          lambda_j = lambda(j);
        end
        if isempty(lambda1)
          lambda1_k = nan;
        else
          lambda1_k = lambda1(k);
        end

        err1 = results(iii).err1;
        err2 = results(iii).err2;
        p1 = results(iii).p1;
        p2 = results(iii).p2;
        cor1 = results(iii).cor1;
        cor2 = results(iii).cor2;
        fprintf('%3d | %6.2f | %6.2f | %10.2f | %10.2f | %10d\n', ...
          icv, lambda_j,lambda1_k,err1,err2,k1);
      end % lam1 loop
    end % lam loop
  end % cv loop
  fprintf('Exit status -- %s (%d iterations)\n', info.message, info.iter);
end % learn_similarity_encoding
