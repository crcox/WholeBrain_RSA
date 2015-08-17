datafile = 's02_avg.mat';
for i = 1:10
  metadata(i).subject = i;
end

params = loadjson('params.json');
datafile = params.data;
[path,fname,ext] = fileparts(datafile);
subjid = sscanf(fname, 's%d');
if ~isempty(subjid)
  subject = find(subjid == [metadata.subject]);
end
metadata = metadata(subject);
