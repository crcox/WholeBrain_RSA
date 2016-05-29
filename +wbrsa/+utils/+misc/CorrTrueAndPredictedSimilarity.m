function [results_out,cor, cor_train] = CorrTrueAndPredictedSimilarity(results)
  n = length(results);

  cor = zeros(1,n);
  cor_train = zeros(1,n);
  for i = 1:n
    r = results(i);
    test = r.cvfilter;
    train = ~test;
    test(r.finalfilter) = [];
    train(r.finalfilter) = [];

    S = r.S(test,test);
    Sz = r.Sz(test,test);
    lt = tril(true(size(S)),-1);
    cor(i) = corr(S(lt),Sz(lt));
    r.cor = cor(i);

    S = r.S(train,train);
    Sz = r.Sz(train,train);
    lt = tril(true(size(S)),-1);
    cor_train(i) = corr(S(lt),Sz(lt));
    r.cor_train = cor_train(i);
    results_out(i) = r;
  end
end
