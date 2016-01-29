function [] = computeGroupStats(results, coordslabel)
  max_xyz = -inf(1,3);
  min_xyz = inf(1,3);
  for i = 1:numel(results)
    z = strcmp({results(i).coords.label}, coordslabel);
    max_xyz = max([max_xyz; results(i).coords.xyz]);
    min_xyz = min([min_xyz; results(i).coords.xyz]);
  end

end
