function WriteAllPermNodestrength(results, outdir, coordslabel)
  if ~exist(outdir, 'dir')
    mkdir(outdir)
  end
  nsubj = numel(results);
  fprintf('Subject: ');
  nchar = 0;
  for i = 1:nsubj
    fprintf(repmat('\b',1,nchar));
    nchar = fprintf('% 4d', i);
    nperm = size(results(i).nodestrength, 1);
    for ii = 1:nperm
      pdir = fullfile(outdir,sprintf('%03d',ii));
      if ~exist(outdir, 'dir');
        mkdir(pdir);
      end
      if results(i).bias
        z = results(i).nz_rows(1:end-1);
      else
        z = results(i).nz_rows;
      end
      val = results(i).nodestrength(ii, z);
      z = strcmp({results(i).coords.label}, coordslabel);
      xyz = results(i).coords(z).xyz;
      z = val > 0;
      val = val(z);
      xyz = xyz(z,:);
      dlmwrite(fullfile(pdir, sprintf('%02d_nodestrength.mni',i)), [xyz,val(:)], ' ');
    end
  end
  fprintf('\n');
end
