function [UzAll,SzAll,nz_rows, p1, p2, train_err, test_err, dist_err] = learn_similarity_encoding(S, V, lambda_try, cvind, holdout, normalize, Gtype, DEBUG, opts)
  [n,d] = size(V);
  V_org = V;

  tt     = 0;
  nt     = n - tt;
  nlam   = length(lambda_try);
  ncv    = max(cvind);
  CVsize = n/ncv;

  if isempty(holdout)
    cvset = 1:ncv;
  else
    cvset = holdout;
  end

  train_err = zeros(nlam,1);
  test_err  = zeros(nlam,1);
  dist_err  = zeros(nlam,1);
  p1        = zeros(nlam,1);
  p2        = zeros(nlam,1);
  nz_rows   = zeros(max(cvind),d,nlam);
  UzAll     = cell(max(cvind),nlam);
  SzAll     = cell(max(cvind),nlam);

  %% preprocessing

  %square root
  [C, r] = sqrt_truncate_r(S, 0.2);

  fprintf('%8s%6s%11s  %11s  %11s  %11s  %11s  %11s  \n', '','lambda','train err','test err', 'p1', 'dist err','p2','n vox')
  for i = cvset
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
          disp('Uz random for debugging...')
          Uz = randn(d,r);
        else
          Uz = Adlas1(V1, C1, lambda,opts);
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
%      lt             = logical(tril(true(nnz(test_set)),0));
%      s              = S(test_set,test_set);
%      sz             = Sz(test_set,test_set);
%      s              = s(lt);
%      sz             = sz(lt);

      p1(i,j)        = trace(corr(S(test_set,:)',Sz(test_set,:)'))/CVsize;
      %p1(i,j)       = corr(s,sz);
      p2(i,j)        = trace(corr(C(test_set,:)',Cz(test_set,:)'))/CVsize;
      train_err(i,j) = norm(S1 - V1*Wz*V1','fro')/norm(S1,'fro');
      test_err(i,j)  = norm(Sz(test_set,:)-S(test_set,:),'fro')/norm(S(test_set,:),'fro');
      dist_err(i,j)  = norm(C(test_set,:) - Cz(test_set,:),'fro')/norm(C(test_set,:),'fro');

      fprintf('%10.2f | %10.2f | %10.2f | %10.2f | %10.2f | %10d\n', ...
        train_err(i,j),test_err(i,j),p1(i,j),dist_err(i,j),p2(i,j),k1);
    end
  end
