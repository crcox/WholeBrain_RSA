function [results,info] = learn_similarity_encoding(S, V, lambda, lambda1, cvind, holdout, normalize, Gtype, DEBUG, options)
  [n,d] = size(V);
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
  nz_rows = zeros(ncv,d,nlam,nlam1);
  UzAll   = cell(ncv,nlam,nlam1);
  SzAll   = cell(ncv,nlam,nlam1);
  if nlam > 1 || nlam1 >1
    nz_rows = squeeze(mat2cell(nz_rows, ncv, d, ones(1,nlam), ones(1,nlam1)));
  end

  %square root
  [C, r] = sqrt_truncate_r(S, 0.2);

  fprintf('%8s%6s%11s  %11s  %11s  %11s  %11s  %11s %11s  \n', '','lambda','test err','train err','p1 test','p1 train','cor test','cor_train','n vox')
  for i = cvset
    CVsize = nnz(cvind==i);
    test_set  = cvind==i;
    train_set = ~test_set;
    fprintf('cv %3d: ', i)

    V = Vorig;
    if normalize == 1
      mm = ones(n,1)*mean(V(train_set,:),1);
      ss = ones(n,1)*std(V(train_set,:),1);
      if all(V(:,end)==1)
        ss(:,d) = ones(n,1);
      end
      V = (V-mm)./ss;
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

        case 'DEBUG'
          Uz = randn(d,r);
          info.message = 'DEBUG'
        end

        UzAll{i,j,k} = Uz;
        fprintf('%6.2f ', lambda)

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
        lt             = logical(tril(true(nnz(test_set)),0));
        s              = S(test_set,test_set);
        sz             = Sz(test_set,test_set);
        st             = St(test_set,test_set);
        s              = s(lt);
        sz             = sz(lt);
        st             = st(lt);

        lt1            = logical(tril(true(nnz(train_set)),0));
        s1             = S(train_set,train_set);
        sz1            = Sz(train_set,train_set);
        st1            = St(train_set,train_set);
        s1             = s1(lt1);
        sz1            = sz1(lt1);
        st1            = st(lt1);

        % Comparison to true S matrix
        p1(i,j,k)       = trace(corr(S(test_set,:)',Sz(test_set,:)'))/nnz(test_set);
        p2(i,j,k)       = trace(corr(S(train_set,:)',Sz(train_set,:)'))/nnz(train_set);
        cor1(i,j,k)     = corr(s,sz); % test
        cor2(i,j,k)     = corr(s1,sz1); % train
        % Comparison to S_trunc = C * C'
        p1t(i,j,k)       = trace(corr(St(test_set,:)',Sz(test_set,:)'))/nnz(test_set);
        p2t(i,j,k)       = trace(corr(St(train_set,:)',Sz(train_set,:)'))/nnz(train_set);
        cor1t(i,j,k)     = corr(st,sz); % test
        cor2t(i,j,k)     = corr(st1,sz1); % train
        % Comparison C
        err1(i,j,k)  = norm(C(test_set,:) - Cz(test_set,:),'fro')/norm(C(test_set,:),'fro');
        err2(i,j,k)  = norm(C(train_set,:) - Cz(train_set,:),'fro')/norm(C(train_set,:),'fro');

        fprintf('%10.2f | %10.2f | %10.2f | %10.2f | %10.2f | %10.2f | %10d\n', ...
          err1(i,j),err2(i,j),p1(i,j),p2(i,j),cor1(i,j),cor2(i,j),k1);

        fprintf('Exit status -- %s (%d iterations)\n', info.message, info.iter);
      end % lam1 loop
    end % lam loop
    results.Uz     = UzAll;
    results.Sz     = SzAll;
    results.nz_rows= nz_rows;
    results.p1     = p1;
    results.p2     = p2;
    results.cor1   = cor1;
    results.cor2   = cor2;
    results.p1t    = p1t;
    results.p2t    = p2t;
    results.cor1t  = cor1t;
    results.cor2t  = cor2t;
    results.err1   = err1;
    results.err2   = err2;
  end % cv loop
end % learn_similarity_encoding
