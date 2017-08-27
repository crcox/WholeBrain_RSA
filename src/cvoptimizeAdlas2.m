function err_growl = optimizeAdlas2(X,Y,c,lamvals)
  X1 = X(test(c),:);
  Y1 = Y(test(c),:);
  S1 = Y1 * Y1';
  X2 = X(training(c),:);
  Y2 = Y(training(c),:);
  lam = lamvals(1);
  lam1 = lamvals(2);
  if all(lamvals > 0) 
    Uz = Adlas2(X2, Y2, lam, lam1, struct());

    Yz1 = X1 * Uz;
    Sz1 = Yz1 * Yz1';
    err_Y = norm(Y1 - Yz1,'fro') / norm(Y1,'fro');
    err_growl = norm(S1 - Sz1,'fro') / norm(S1,'fro');
  else
    err_growl = inf;
  end
end
