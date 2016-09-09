function results_tbl = mystruct2tbl(results)
  fields = fieldnames(results);
  scalarfields = false(1,numel(fieldnames(results)));
  for i = 1:numel(fields)
    f = fields{i};
    x = results(1).(f);
    if (isscalar(x) || ischar(x)) && ~isstruct(x)
      scalarfields(i) = 1;
    end
  end
  z = cellfun(@isempty, {results.subject});
  results_tbl = struct2table(rmfield(results(~z),fields(~scalarfields)));
end