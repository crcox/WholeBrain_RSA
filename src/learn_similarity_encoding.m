function [results,info] = learn_similarity_encoding(S, V, Gtype, varargin)
  p = inputParser();
  addRequired(p,'S');
  addRequired(p,'V');
  addRequired(p,'Gtype');
  addParameter(p , 'tau'       , 0.2);
  addParameter(p , 'lambda'    , []);
  addParameter(p , 'lambda1'   , []);
  addParameter(p , 'cvind'     , []);
  addParameter(p , 'cvholdout' , []);
  addParameter(p , 'normalize' , []);
  addParameter(p , 'LambdaSeq' , []);
  addParameter(p , 'DEBUG'     , false);
  addParameter(p , 'AdlasOpts' , struct());
  addParameter(p , 'SmallFootprint' ,false);
  parse(p, S, V, Gtype, varargin{:});

  S         = p.Results.S;
  V         = p.Results.V;
  Gtype     = p.Results.Gtype;
  tau       = p.Results.tau;
  lambda    = p.Results.lambda;
  lambda1   = p.Results.lambda1;
  cvind     = p.Results.cvind;
  holdout   = p.Results.cvholdout;
  normalize = p.Results.normalize;
  LambdaSeq = p.Results.LambdaSeq;
  DEBUG     = p.Results.DEBUG;
  options   = p.Results.AdlasOpts;
  SMALL     = p.Results.SmallFootprint;

  if strcmp(Gtype, {'grOWL','grOWL2'});
    assert(~isempty(LambdaSeq),'A LambdaSeq type (linear or exponential) must be set when using grOWL*');
  end

  [n,d] = size(V);
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

  %ncv = max(cvind);
  if isempty(holdout)
    cvset = 1:ncv;
  else
    cvset = holdout;
  end
  ncv = numel(cvset);

  % Define results structure
  results.Uz = []
  results.Cz = []
  results.Sz = []
  results.nz_rows =  [];
  results.subject =  [];
  results.cvholdout = [];
  results.finalholdout = [];
  results.lambda = [];
  results.lambda1 = [];
  results.LambdaSeq = [];
  results.Gtype = [];
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

  fprintf('%8s%6s%11s %11s  %11s  %11s  %11s  %11s  %11s %11s  \n', '','lam','lam1','test err','train err','p1 test','p1 train','cor test','cor train','n vox')

  iii = 0; % index into 1-D results structure.
  for i = 1:ncv
    icv = cvset(i);
    CVsize = nnz(cvind==icv);
    test_set  = cvind==icv;
    train_set = ~test_set;
    fprintf('cv %3d: ', icv)

    V = Vorig;

    % remove bias unit for normalization
    if BIASUNIT
      V(:,end) = [];
    end

    % normalize
    switch normalize
    case 'zscore_train'
      mm = mean(V(train_set,:),1);
      ss = std(V(train_set,:),0,1);
    case 'zscore'
      mm = mean(V,1);
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
    V = bsxfun(@rdivide,V, ss);

    % Reinsert bias term
    if BIASUNIT
      V(:,end+1) = 1;
    end

    C1 = C(train_set,:);
    V1 = V(train_set,:);

    S1 = S(train_set, train_set);
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
          switch Gtype
          case 'L1L2'
            lam1 = lambda1(j);
            if all(lam==0)
              Uz = pinv(V1)*C1;
            else
              [Uz, info] = Adlas1(V1, C1, lam1, options);
            end

          case 'grOWL'
            switch LambdaSeq
            case 'linear'
              lam = lam(j)*(d:-1:1)/d;
            case 'exponential'
              lam = lam(j)*sqrt(2*log((d*ones(1,d))./(1:d)));
            end
            if all(lam==0)
              Uz = pinv(V1)*C1;
            else
              [Uz, info] = Adlas1(V1, C1, lam, options);
            end

          case 'grOWL2'
            lam1 = lambda1(k);
            switch LambdaSeq
            case 'linear'
              lam = lam(j)*(d:-1:1)/d;
            case 'exponential'
              lam = lam(j)*sqrt(2*log((d*ones(1,d))./(1:d)));
            case 'inf' % needs both lam and lam1
              lam = [lam+lam1, repmat(lam,1,d-1)];
            end
            if all(lam==0)
              Uz = pinv(V1)*C1;
            else
              [Uz, info] = Adlas1(V1, C1, lam, options);
            end
          end
        end

        % KEY:
        % Our models assume that S = V*W*V'
        % V  : The fMRI data.
        % W  : A matrix of weights.
        % S  : The true similarity matrix.
        % Sz : The predicted similarity matrix.
        % C  : Essentially the first r principle components of S.
        % Cz : The predicted C.
        % St : The approximated S, reconstructed from actual C.
        % Uz : The estimated voxel weights, with a r weights per voxel.

        % Prepare to evaluate solutions
        k1 = nnz(any(Uz,2));
        Wz = Uz*Uz';
        Sz = V*Wz*V';
        Cz = V*Uz;
        St = C*C';

        lt  = logical(tril(true(nnz(test_set)),0));
        s   = S(test_set,test_set);
        sz  = Sz(test_set,test_set);
        st  = St(test_set,test_set);
        s   = s(lt);
        sz  = sz(lt);
        st  = st(lt);

        lt2 = logical(tril(true(nnz(train_set)),0));
        s2  = S(train_set,train_set);
        sz2 = Sz(train_set,train_set);
        st2 = St(train_set,train_set);
        s2  = s1(lt1);
        sz2 = sz1(lt1);
        st2 = st1(lt1);

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
        results(iii).Gtype = Gtype;
        results(iii).bias = BIASUNIT;
        results(iii).normalize = normalize;
        results(iii).nzv = k1;
        % Comparison to true S matrix
        results(iii).p1      = trace(corr(S(test_set,:)',Sz(test_set,:)'))/nnz(test_set);
        results(iii).p2      = trace(corr(S(train_set,:)',Sz(train_set,:)'))/nnz(train_set);
        results(iii).cor1    = corr(s,sz); % test
        results(iii).cor2    = corr(s2,sz2); % train
        % Comparison to S_trunc = C * C'
        results(iii).p1t     = trace(corr(St(test_set,:)',Sz(test_set,:)'))/nnz(test_set);
        results(iii).p2t     = trace(corr(St(train_set,:)',Sz(train_set,:)'))/nnz(train_set);
        results(iii).cor1t   = corr(st,sz); % test
        results(iii).cor2t   = corr(st2,sz2); % train
        % Comparison to C
        results(iii).err1    = norm(C(test_set,:) - Cz(test_set,:),'fro')/norm(C(test_set,:),'fro');
        results(iii).err2    = norm(C(train_set,:) - Cz(train_set,:),'fro')/norm(C(train_set,:),'fro');
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

        fprintf('%6.2f | %6.2f | %10.2f | %10.2f | %10.2f | %10.2f | %10.2f | %10.2f | %10d\n', ...
          lambda_j,lambda1_k,err1(i,j,k),err2(i,j,k),p1(i,j,k),p2(i,j,k),cor1(i,j,k),cor2(i,j,k),k1);

        fprintf('Exit status -- %s (%d iterations)\n', info.message, info.iter);
      end % lam1 loop
    end % lam loop
  end % cv loop
end % learn_similarity_encoding
