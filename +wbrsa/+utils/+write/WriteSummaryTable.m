function WriteSummaryTable(results,params)
  ModelTypes = {'use_random_cor','use_random_inner','use_shuffled_cor','use_shuffled_inner','real'};
  DataTypes = {'userandom','useshuffled','real'};
  Normalize = [true,false];
  Bias = [true,false];

  fid = fopen('summary.csv', 'w');
  for i = 1:length(ModelTypes)
    mt = ModelTypes{i};
    for j = 1:length(DataTypes)
      dt = DataTypes{j};
      for k = 1:length(Normalize)
        norm = Normalize(k);
        for m = 1:length(Bias)
          bias = Bias(m);
          z = strcmp(mt,{params.SanityCheckModel}) & ...
              strcmp(dt,{params.SanityCheckData}) & ...
              [params.normalize] == norm & ...
              [params.bias] == bias;
          n = max([results(z).cvind])-1;
          x1 = mean([results(z).p1]);
          s = std([results(z).p1])/sqrt(n);
          c1 = s*2;

          if isfield(results,'cor')
            x2 = mean([results(z).cor]);
            s = std([results(z).cor])/sqrt(n);
            c2 = s*2;
            fprintf(fid,'%s,%s,%d,%d,%6.3f,%6.3f,%6.3f,%6.3f\n',mt,dt,int32(norm),int32(bias),x1,c1,x2,c2);
          else
            fprintf(fid,'%s,%s,%d,%d,%6.3f,%6.3f\n',mt,dt,int32(norm),int32(bias),x1,c1);
          end

        end
      end
    end
  end
  fclose(fid);
end
