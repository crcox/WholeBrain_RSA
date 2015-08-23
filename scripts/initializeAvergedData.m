% INITIALIZE AVERAGED DATA 
% Note that in this scheme, the data are replicated to match the similarity
% structures.

% pathtorepo = '~/src/WholeBrain_RSA/'; % lab machine
pathtorepo = 'C:/Users/chris/Documents/WholeBrain_RSA/'; % home machine
addpath(fullfile(pathtorepo,'src'),fullfile(pathtorepo,'util'))
% datadir_in = '~/data/Manchester/data/mat/fromRick';
% datadir_out = '~/data/Manchester/WholeBrain_RSA/data/avg';
datadir_in = './UrvashiData';
datadir_out = './DataTest';

%% Define metadata
% number of subjects
n = 23;

% presentation modality
% There are 37 visual and 37 auditory stimuli that are paired. So there are
% 74 unique stimuli, that relate to 37 concepts. Stimuli are ordered:
% 37 visual , first presentation
% 37 audio  , first presentation
% 37 visual , second presentation
% 37 audio  , second presentation
% etc ...
% Each stimulus is presented 4 times.
stimcode = repmat(1:74, 1, 4);
conceptcode = repmat(1:37, 1, 8);
visual = stimcode < 38;
audio = stimcode > 37;

%% Compute outliers and averages
data_varname_in = 'X';
filters = struct();
fmt_header = '% 9s% 6s% 8s% 10s% 6s\n';
fmt_err = '% 9d% 6s% 8s% 10s% 6s\n';
fmt = '% 9d% 6s% 8s% 10s% 6s';

fprintf(fmt_header, 'subject', 'load', 'filter', 'outliers', 'save');
for i = 1:n;
  loaded = '';
  outliers_identified = '';
  filtered = '';
  saved = '';
  strlen = fprintf(fmt, i, loaded, outliers_identified, filtered, saved);
  
  name_in = sprintf('subject%02d.mat',i);
  name_out = sprintf('s%02d_avg.mat',i);
  datapath_in = fullfile(datadir_in, name_in);
  datapath_out = fullfile(datadir_out, name_out);
  
  % Load
  try
    StagingContainer = load(datapath_in, data_varname_in);
    Data.raw = StagingContainer.(data_varname_in);
  catch ME
    loaded = 'ERROR';
    fprintf(repmat('\b',1,strlen));
    fprintf(fmt_err, i, loaded, outliers_identified, filtered, saved); 
    rethrow(ME);
  end
  loaded = 'ok';
  fprintf(repmat('\b',1,strlen));
  strlen = fprintf(fmt, i, loaded, outliers_identified, filtered, saved);
  clear StagingContainer;
  
  % Filter and Aggregate
  try
    Averaged.visual = averageRepeatedTrials(Data.raw(visual,:), stimcode(visual));
    Averaged.audio = averageRepeatedTrials(Data.raw(audio,:), stimcode(audio));
    Averaged.semantic = averageRepeatedTrials(Data.raw, conceptcode);
  catch ME
    filtered = 'ERROR';
    fprintf(repmat('\b',1,strlen));
    fprintf(fmt_err, i, loaded, filtered, outliers_identified, saved);
    rethrow(ME);
  end
  filtered = 'ok';
  fprintf(repmat('\b',1,strlen));
  strlen = fprintf(fmt, i, loaded, filtered, outliers_identified, saved);
  clear Data;
  
  % Identify outliers
  try
    [~,filters(i).visual] = removeOutliers(Averaged.visual);
    [~,filters(i).audio] = removeOutliers(Averaged.audio);
    [~,filters(i).semantic] = removeOutliers(Averaged.semantic);
  catch ME
    outliers_identified = 'ERROR';
    fprintf(repmat('\b',1,strlen));
    fprintf(fmt_err, i, loaded, filtered, outliers_identified, saved);
    rethrow(ME);
  end
  outliers_identified = 'ok';
  fprintf(repmat('\b',1,strlen));
  strlen = fprintf(fmt, i, loaded, filtered, outliers_identified, saved);

  % save
  try
    save(datapath_out, '-struct', 'Averaged'); 
  catch ME
    saved = 'ERROR';
    fprintf(repmat('\b',1,strlen));
    fprintf(fmt_err, i, loaded, filtered, outliers_identified, saved);
    rethrow(ME);
  end
  saved = 'ok';
  fprintf(repmat('\b',1,strlen));
  strlen = fprintf(fmt, i, loaded, filtered, outliers_identified, saved);
  fprintf('\n');
end

%% Compose metadata
metadata = struct();
for i = 1:n
  metadata(i).datacode = 'avg';
  metadata(i).subject = i;
  metadata(i).visual = visual;
  metadata(i).audio = audio;
  metadata(i).conceptcode = conceptcode;
  metadata(i).stimcode = stimcode;
  %
  metadata(i).filter(1).subset = 'visual';
  metadata(i).filter(1).label = 'rowfilter';
  metadata(i).filter(1).dimension = 1;
  metadata(i).filter(1).filter = filters(i).visual.words;
  metadata(i).filter(1).notes = 'Subset by modality, average over stimulus repetitions.';
  %
  metadata(i).filter(2).subset = 'visual';
  metadata(i).filter(2).label = 'colfilter';
  metadata(i).filter(2).dimension = 2;
  metadata(i).filter(2).filter = filters(i).visual.voxels;
  metadata(i).filter(2).notes = 'Subset by modality, average over stimulus repetitions.';
  %
  metadata(i).filter(3).subset = 'audio';
  metadata(i).filter(3).label = 'rowfilter';
  metadata(i).filter(3).dimension = 1;
  metadata(i).filter(3).filter = filters(i).audio.words;
  metadata(i).filter(3).notes = 'Subset by modality, average over stimulus repetitions.';
  %
  metadata(i).filter(4).subset = 'audio';
  metadata(i).filter(4).label = 'colfilter';
  metadata(i).filter(4).dimension = 2;
  metadata(i).filter(4).filter = filters(i).audio.voxels;
  metadata(i).filter(4).notes = 'Subset by modality, average over stimulus repetitions.';
  %
  metadata(i).filter(5).subset = 'semantic';
  metadata(i).filter(5).label = 'rowfilter';
  metadata(i).filter(5).dimension = 1;
  metadata(i).filter(5).filter = filters(i).semantic.words;
  metadata(i).filter(5).notes = 'Average over stimulus repetitions and modalities.';
  %
  metadata(i).filter(6).subset = 'semantic';
  metadata(i).filter(6).label = 'colfilter';
  metadata(i).filter(6).dimension = 2;
  metadata(i).filter(6).filter = filters(i).semantic.voxels;
  metadata(i).filter(6).notes = 'Average over stimulus repetitions and modalities.';

end
metapath_out = fullfile(datadir_out, 'metadata_avg.mat');
save(metapath_out, 'metadata');

%% Copy similarity structures into place
% They are already there, so I am not going to over engineer this...

%% Copy coordinates into place
% They are already there, so I am not going to over engineer this...

%% Pre-generate CV schemes
% There are 37 unique concepts.
CV = defineCVBlocks(1:37, 'folds', 9, 'schemes', 11);
CV(:,12) = defineCVBlocks(1:37, 'LOO', true);
save(fullfile(datadir_out, 'CV_schemes_avg.mat'), 'CV');