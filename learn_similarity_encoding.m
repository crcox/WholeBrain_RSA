function [nz_rows, p1, p2, train_err, test_err, dist_err] = learn_similarity_encoding(S, V, lambda_try, cvind, normalize, Gtype, DEBUG)
  [n,d] = size(V);
  V_org = [V, ones(n,1)];
  d = d + 1;

  tt = 0;
  nt = n - tt;
  nlam = length(lambda_try);
  ncv = max(cvind);
  CVsize = n/ncv;

  train_err = zeros(nlam,1);
  test_err = zeros(nlam,1);
  dist_err = zeros(nlam,1);
  p1 = zeros(nlam,1);
  p2 = zeros(nlam,1);
  nz_rows = zeros(max(cvind),d,nlam);
  UzAll = cell(max(cvind),nlam);

  %% preprocessing

  %square root
  % Higher values of the second parameter will lead to lower-rank matrices
  % because it permits more error between the truncated and the original
  % matrix.
  [C, r] = sqrt_truncate_r(S, 0.2);

  for i = 1:ncv
    V = V_org;
    test_set  = cvind==i;
    train_set = ~test_set;
    fprintf('cv% 3d: ', i)

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
        Uz = Adlas1(V1, C1, lambda);
      else
        lambda =  lambda_try(j);
        if DEBUG
          Uz = randn(d,r);
        else
          Uz = Adlas1(V1, C1, lambda);
        end
      end
      UzCell{i,j} = Uz;
      fprintf('%.2f ', lambda)

      k1  = sum(any(Uz,2));
      Wz = Uz*Uz';
      Sz = V*Wz*V';
      Cz = V*Uz;

      %store results
      nz_rows(i,:,j) = any(Uz,2);
      p1(i,j) = trace(corr(S(test_set,:)',Sz(test_set,:)'))/CVsize;
      p2(i,j) = trace(corr(C(test_set,:)',Cz(test_set,:)'))/CVsize;
      train_err(i,j) = norm(S1 - V1*Wz*V1','fro')/norm(S1,'fro');
      test_err(i,j)  = norm(Sz(test_set,:)-S(test_set,:),'fro')/norm(S(test_set,:),'fro');
      dist_err(i,j)  = norm(C(test_set,:) - Cz(test_set,:),'fro')/norm(C(test_set,:),'fro');

      fprintf('% 3d | % 10.2f | % 10.2f | % 10.2f | % 10.2f | % 10.2f | % 10.2f | % 10.2f\n',i,lambda_try(j),train_err(i,j),test_err(i,j),p1(i,j),dist_err(i,j),p2(i,j),k1);
    end
  end
