function Y = selectTargets(metadata, type, label, source, metric, rowfilter)
  Y = cell(1, numel(metadata));
  if ~iscell(rowfilter);
    rowfilter = {rowfilter};
  end
  for i = 1:numel(metadata);
    T = selectbyfield(metadata(i).targets, 'label', label, 'type', type, 'sim_source', source, 'sim_metric', metric);
    if numel(T) == 1
      switch lower(type)
        case 'similarity'
          Y{i} = T.target(rowfilter{i}, rowfilter{i});
        case 'embedding'
          Y{i} = T.target(rowfilter{i}, :);
        case 'category'
          Y{i} = T.target(rowfilter{i});
      end
    elseif numel(T) == 0;
        error('crcox:NoTargetsSelected', 'There are no targets in the metadata structure that satisfy the criteria.');
    else
        error('crcox:ManyTargetsSelected', 'There are many targets in the metadata structure that satisfy the criteria.');
    end
  end
end
