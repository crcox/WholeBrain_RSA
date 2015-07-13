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

  ncv = max(cvind);
  if isempty(holdout)
    cvset = 1:ncv;
  else
    cvset = holdout;
  end

  err1    = zeros(ncv,nlam,nlam1);
  err2    = zeros(ncv,nlam,nlam1);
  cor1    = zeros(ncv,nlam,nlam1);
  cor2    = zeros(ncv,nlam,nlam1);
  p1      = zeros(ncv,nlam,nlam1);
  p2      = zeros(ncv,nlam,nlam1);
  cor1t   = zeros(ncv,nlam,nlam1);
  cor2t   = zeros(ncv,nlam,nlam1);
  p1t     = zeros(ncv,nlam,nlam1);
  p2t     = zeros(ncv,nlam,nlam1);
  nz_rows = false(ncv,d,nlam,nlam1);
  UzAll   = cell(ncv,nlam,nlam1);
  SzAll   = cell(ncv,nlam,nlam1);
  if nlam > 1 || nlam1 > 1
    nz_rows = squeeze(mat2cell(nz_rows, ncv, d, ones(1,nlam), ones(1,nlam1)));
  end

  %square root
  [C, r] = sqrt_truncate_r(S, tau);

  fprintf('%8s%6s%11s %11s  %11s  %11s  %11s  %11s  %11s %11s  \n', '','lam','lam1','test err','train err','p1 test','p1 train','cor test','cor train','n vox')
  for i = cvset
    CVsize = nnz(cvind==i);
    test_set  = cvind==i;
    train_set = ~test_set;
    fprintf('cv %3d: ', i)

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
    case '2norm'
      mm = mean(V,1);
      ss = norm(V);
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
            [Uz,info] = Adlas1(V1, C1, lam1, options);

          case 'grOWL'
            switch LambdaSeq
            case 'linear'
              lam = lam(j)*(d:-1:1)/d;
            case 'exponential'
              lam = lam(j)*sqrt(2*log((d*ones(1,d))./(1:d)));
            end
            [Uz,info] = Adlas1(V1, C1, lam, options);

          case 'grOWL2'
            switch LambdaSeq
            case 'linear'
              lam = lam(j)*(d:-1:1)/d;
            case 'exponential'
              lam = lam(j)*sqrt(2*log((d*ones(1,d))./(1:d)));
            end
            lam1 = lambda1(k);
            [Uz, info] = Adlas2(V1, C1, lam, lam1, options);
          end
        end

        UzAll{i,j,k} = Uz;

        k1 = nnz(any(Uz,2));
        Wz = Uz*Uz';
        Sz = V*Wz*V';
        Cz = V*Uz;
        St = C*C';

        SzAll{i,j,k} = Sz;

        %store results
        if nlam > 1 || nlam1 >1
          nz_rows{j,k}(i,:) = any(Uz,2);
        else
          nz_rows(i,:) = any(Uz,2);
        end
        lt  = logical(tril(true(nnz(test_set)),0));
        s   = S(test_set,test_set);
        sz  = Sz(test_set,test_set);
        st  = St(test_set,test_set);
        s   = s(lt);
        sz  = sz(lt);
        st  = st(lt);

        lt1 = logical(tril(true(nnz(train_set)),0));
        s1  = S(train_set,train_set);
        sz1 = Sz(train_set,train_set);
        st1 = St(train_set,train_set);
        s1  = s1(lt1);
        sz1 = sz1(lt1);
        st1 = st1(lt1);

        % Comparison to true S matrix
        p1(i,j,k)    = trace(corr(S(test_set,:)',Sz(test_set,:)'))/nnz(test_set);
        p2(i,j,k)    = trace(corr(S(train_set,:)',Sz(train_set,:)'))/nnz(train_set);
        cor1(i,j,k)  = corr(s,sz); % test
        cor2(i,j,k)  = corr(s1,sz1); % train
        % Comparison to S_trunc = C * C'
        p1t(i,j,k)   = trace(corr(St(test_set,:)',Sz(test_set,:)'))/nnz(test_set);
        p2t(i,j,k)   = trace(corr(St(train_set,:)',Sz(train_set,:)'))/nnz(train_set);
        cor1t(i,j,k) = corr(st,sz); % test
        cor2t(i,j,k) = corr(st1,sz1); % train
        % Comparison to C
        err1(i,j,k)  = norm(C(test_set,:) - Cz(test_set,:),'fro')/norm(C(test_set,:),'fro');
        err2(i,j,k)  = norm(C(train_set,:) - Cz(train_set,:),'fro')/norm(C(train_set,:),'fro');

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
  if ~SMALL
    results.Uz = UzAll;
    results.Sz = SzAll;
  end
  if ~iscell(nz_rows);
    results.nz_rows = nz_rows(cvset,:,:);
  else
    results.nz_rows = nz_rows;
  end
  results.p1      = p1(cvset,:,:);
  results.p2      = p2(cvset,:,:);
  results.cor1    = cor1(cvset,:,:);
  results.cor2    = cor2(cvset,:,:);
  results.p1t     = p1t(cvset,:,:);
  results.p2t     = p2t(cvset,:,:);
  results.cor1t   = cor1t(cvset,:,:);
  results.cor2t   = cor2t(cvset,:,:);
  results.err1    = err1(cvset,:,:);
  results.err2    = err2(cvset,:,:);
  results.iter    = info.iter;
end % learn_similarity_encoding
