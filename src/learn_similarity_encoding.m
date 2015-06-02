function [results,info] = learn_similarity_encoding(S, V, lambda_try, cvind, holdout, normalize, Gtype, DEBUG, options)
  [n,d] = size(V);
  Vorig = V;

  tt     = 0;
  nt     = n - tt;
  nlam   = length(lambda_try);
  ncv    = max(cvind);

  if isempty(holdout)
    cvset = 1:ncv;
  else
    cvset = holdout;
  end

  err1    = zeros(nlam,1);
  err2    = zeros(nlam,1);
  cor1    = zeros(nlam,1);
  cor2    = zeros(nlam,1);
  p1      = zeros(nlam,1);
  p2      = zeros(nlam,1);
  cor1t    = zeros(nlam,1);
  cor2t    = zeros(nlam,1);
  p1t      = zeros(nlam,1);
  p2t      = zeros(nlam,1);
  nz_rows = zeros(max(cvind),d,nlam);
  UzAll   = cell(max(cvind),nlam);
  SzAll   = cell(max(cvind),nlam);

  %% preprocessing

  %square root
  [C, r] = sqrt_truncate_r(S, 0.2);

  fprintf('%8s%6s%11s  %11s  %11s  %11s  %11s  %11s %11s  \n', '','lambda','test err','train err', 'p1 test', 'p1 train','cor test','cor_train','n vox')
  for i = cvset
    CVsize = nnz(cvind==i);
    V = Vorig;
    test_set  = cvind==i;
    train_set = ~test_set;
    fprintf('cv %3d: ', i)

    if normalize == 1
      mm = ones(n,1)*mean(Vorig(train_set,1);
      ss = ones(n,1)*std(Vorig(train_set,1);
      if all(V(:,end)==1)
        ss(:,d) = ones(n,1);
      end
      V = (Vorig-mm)./ss;
    end

    C1 = C(train_set,:);
    V1 = V(train_set,:);

    S1 = S(train_set, train_set);
    for j = 1:length(lambda_try)

      switch Gtype
      case 'grOWL'
        lambda = lambda_try(j)*(d:-1:1)/d;
        [Uz,info] = Adlas1(V1, C1, lambda, options);

      case 'grOWL2'
        lambda = lambda_try(j)*(d:-1:1)/d;
        lambda1 = lambda_try(j);
        [Uz, info] = Adlas2(V1, C1, lambda, lambda1, options);

      case 'L1L2'
        lambda = lambda_try(j);
        [Uz,info] = Adlas1(V1, C1, lambda, options);

      case 'DEBUG'
        Uz = randn(d,r);
        info.message = 'DEBUG'
      end

      UzAll{i,j} = Uz;
      fprintf('%6.2f ', lambda)

      k1 = nnz(any(Uz,2));
      Wz = Uz*Uz';
      Sz = V*Wz*V';
      Cz = V*Uz;
      St = C*C';

      SzAll{i,j} = Sz;

      %store results
      nz_rows(i,:,j) = any(Uz,2);
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
      p1(i,j)       = trace(corr(S(test_set,:)',Sz(test_set,:)'))/nnz(test_set);
      p2(i,j)       = trace(corr(S(train_set,:)',Sz(train_set,:)'))/nnz(train_set);
      cor1(i,j)     = corr(s,sz); % test
      cor2(i,j)     = corr(s1,sz1); % train
      % Comparison to S_trunc = C * C'
      p1t(i,j)       = trace(corr(St(test_set,:)',Sz(test_set,:)'))/nnz(test_set);
      p2t(i,j)       = trace(corr(St(train_set,:)',Sz(train_set,:)'))/nnz(train_set);
      cor1t(i,j)     = corr(st,sz); % test
      cor2t(i,j)     = corr(st1,sz1); % train
      % Comparison C
      err1(i,j)  = norm(C(test_set,:) - Cz(test_set,:),'fro')/norm(C(test_set,:),'fro');
      err2(i,j)  = norm(C(train_set,:) - Cz(train_set,:),'fro')/norm(C(train_set,:),'fro');

      fprintf('%10.2f | %10.2f | %10.2f | %10.2f | %10.2f | %10.2f | %10d\n', ...
        err1(i,j),err2(i,j),p1(i,j),p2(i,j),cor1(i,j),cor2(i,j),k1);

      fprintf('Exit status -- %s (%d iterations)\n', info.message, info.iter);
    end
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
  end
