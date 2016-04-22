function err_L1L2 = optimizeAdlas1(X,Y,c,lam1)
  X1 = X(test(c),:);
  Y1 = Y(test(c),:);
  S1 = Y1 * Y1';
  X2 = X(training(c),:);
  Y2 = Y(training(c),:);

  Uz = Adlas1(X2, Y2, lam1);

  Yz1 = X1 * Uz;
  Sz1 = Yz1 * Yz1';
  err_Y = norm(Y1 - Yz1,'fro')/norm(Y1,'fro');
  err_L1L2 = norm(S1 - Sz1,'fro')/norm(S1,'fro');
end
  
