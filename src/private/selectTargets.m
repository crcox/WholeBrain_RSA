function Y = selectTargets(metadata, type, label, source, metric, rowfilter)
  Y = cell(1, numel(metadata));
  if numel(metadata) == 1;
    rowfilter = {rowfilter};
  end
  for i = 1:numel(metadata);
    t = metadata(i).targets;
    z = all([strcmpi(type,{t.type});strcmpi(label,{t.label});strcmpi(source,{t.sim_source});strcmpi(metric,{t.sim_metric})]);
    if strcmpi(type, 'similarity')
      Y{i} = metadata(i).targets(z).target(rowfilter{i}, rowfilter{i});
    else
      Y{i} = metadata(i).targets(z).target(rowfilter{i});
    end
  end
end
