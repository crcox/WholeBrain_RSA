function [UzAll,SzAll,nz_rows, p1, p2, cor1, cor2, err1, err2 ] = learn_similarity_encoding(S, V, lambda_try, cvind, holdout, normalize, Gtype, DEBUG, opts)
  [n,d] = size(V);
  V_org = V;

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
  nz_rows = zeros(max(cvind),d,nlam);
  UzAll   = cell(max(cvind),nlam);
  SzAll   = cell(max(cvind),nlam);

  %% preprocessing

  %square root
  [C, r] = sqrt_truncate_r(S, 0.2);

  fprintf('%8s%6s%11s  %11s  %11s  %11s  %11s  %11s %11s  \n', '','lambda','test err','train err', 'p1 test', 'p1 train','cor test','cor_train','n vox')
  for i = cvset
    CVsize = nnz(cvind==i);
    V = V_org;
    test_set  = cvind==i;
    train_set = ~test_set;
    fprintf('cv %3d: ', i)

    % CRC: Normalize only based on the training set.
    if normalize == 1
      mm = zeros(n,d);
      ss = ones(n,1)*std(V,1);
      ss(:,d) = ones(n,1);
      V = (V-mm)./ss;
    end

    C1 = C(train_set,:);
    V1 = V(train_set,:);

    S1 = S(train_set, train_set);
    for j = 1:length(lambda_try)
      if strcmp(Gtype, 'grOWL')
        lambda = lambda_try(j)*(d:-1:1)/d;
        Uz = Adlas1(V1, C1, lambda,opts);
      else
        lambda =  lambda_try(j);
        if DEBUG
          Uz = randn(d,r);
          info.message = 'DEBUG'
        else
          [Uz,info] = Adlas1(V1, C1, lambda,opts);
        end
      end
      UzAll{i,j} = Uz;
      fprintf('%6.2f ', lambda)

      k1 = nnz(any(Uz,2));
      Wz = Uz*Uz';
      Sz = V*Wz*V';
      Cz = V*Uz;

      SzAll{i,j} = Sz;
%      S_nd = reshape(S(~eye(n)),n,n-1);
%      Sz_nd = reshape(Sz(~eye(n)),n,n-1);

      %store results
      nz_rows(i,:,j) = any(Uz,2);
      lt             = logical(tril(true(nnz(test_set)),0));
      s              = S(test_set,test_set);
      sz             = Sz(test_set,test_set);
      s              = s(lt);
      sz             = sz(lt);

      lt1            = logical(tril(true(nnz(train_set)),0));
      s1             = S(train_set,train_set);
      sz1            = Sz(train_set,train_set);
      s1             = s1(lt1);
      sz1            = sz1(lt1);

      p1(i,j)       = trace(corr(S(test_set,:)',Sz(test_set,:)'))/nnz(test_set);
      p2(i,j)       = trace(corr(S(train_set,:)',Sz(train_set,:)'))/nnz(train_set);
      cor1(i,j)     = corr(s,sz); % test
      cor2(i,j)     = corr(s1,sz1); % train
%      train_err(i,j) = norm(S1 - V1*Wz*V1','fro')/norm(S1,'fro');
%      test_err(i,j)  = norm(Sz(test_set,:)-S(test_set,:),'fro')/norm(S(test_set,:),'fro');
      err1(i,j)  = norm(C(test_set,:) - Cz(test_set,:),'fro')/norm(C(test_set,:),'fro');
      err2(i,j)  = norm(C(train_set,:) - Cz(train_set,:),'fro')/norm(C(train_set,:),'fro');

      fprintf('%10.2f | %10.2f | %10.2f | %10.2f | %10.2f | %10.2f | %10d\n', ...
        err1(i,j),err2(i,j),p1(i,j),p2(i,j),cor1(i,j),cor2(i,j),k1);

      fprintf('Exit status -- %s\n', info.message);
    end
  end
