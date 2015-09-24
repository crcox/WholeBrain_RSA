addpath('~/src/WholeBrain_RSA/src');

% Based of the avg directory
srcdir = 'avg';
destdir = 'avgVisAud';

% Load metadata from srcdir
tmp = load(fullfile(srcdir, sprintf('metadata_%s.mat', srcdir)), 'metadata');
metadata = tmp.metadata;
clear tmp;

% Update high-level metdata structure
[metadata.datacode] = deal(destdir);

% Load and concatenate subject data. Update filters.
for i = 1:23
  filename_in = sprintf('s%02d_%s.mat', i, srcdir);
  filepath_in = fullfile(srcdir, filename_in);
  filename_out = sprintf('s%02d_%s.mat', i, destdir);
  filepath_out = fullfile(destdir, filename_out);
  load(filepath_in, 'visual', 'audio');
  VisAud = [visual;audio];
  [~,filter] = removeOutliers(VisAud);
  metadata(i).filter(1).subset = 'VisAud';
  metadata(i).filter(1).label = 'rowfilter';
  metadata(i).filter(1).filter = filter.words;
  metadata(i).filter(1).notes = 'Concatenate averaged visual and audio trials.';
  metadata(i).filter(2).subset = 'VisAud';
  metadata(i).filter(2).label = 'colfilter';
  metadata(i).filter(2).filter = filter.voxels;
  metadata(i).filter(2).notes = 'Concatenate averaged visual and audio trials.';
  metadata(i).filter(3:end) = [];
  save(filepath_out,'VisAud');
end

% Save metadata
save(fullfile(destdir, sprintf('metadata_%s.mat', destdir)), 'metadata');

% Load and concatenate similarity matrices
filename_in = sprintf('model_semantic_%s.mat', srcdir);
filepath_in = fullfile(srcdir, filename_in);
filename_out = sprintf('model_semantic_%s.mat', destdir);
filepath_out = fullfile(destdir, filename_out);
load(filepath_in,'cosine');
cosine = repmat(cosine,2,2);
save(filepath_out, 'cosine');

filename_in = sprintf('model_visual_%s.mat', srcdir);
filepath_in = fullfile(srcdir, filename_in);
filename_out = sprintf('model_visual_%s.mat', destdir);
filepath_out = fullfile(destdir, filename_out);
load(filepath_in,'earthmover');
earthmover = repmat(earthmover,2,2);
save(filepath_out, 'earthmover');

% Load and concatenate CV schemes
filename_in = sprintf('CV_schemes_%s.mat', srcdir);
filepath_in = fullfile(srcdir, filename_in);
filename_out = sprintf('CV_schemes_%s.mat', destdir);
filepath_out = fullfile(destdir, filename_out);
load(filepath_in,'CV');
CV = repmat(CV,2,1);
save(filepath_out, 'CV');

% copy over coords
copyfile(fullfile(srcdir,'coords.mat'), fullfile(destdir,'coords.mat'));
