% INITIALIZE REPLICATED DATA
% Note that in this scheme, the similarity structures are replicated to
% match the data.

%pathtorepo = '~/src/WholeBrain_RSA/'; % lab machine
pathtorepo = 'C:/Users/chris/Documents/WholeBrain_RSA/'; % home machine
addpath(fullfile(pathtorepo,'src'),fullfile(pathtorepo,'util'))
% datadir_in = '~/data/Manchester/data/mat/fromRick';
% datadir_out = '~/data/Manchester/WholeBrain_RSA/data/avg';
datadir_in = './UrvashiData';
datadir_out = './DataTest';

%% Define metadata
% number of subjects
n = 1;

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

%% Compute outliers
data_varname_in = 'X';
filters = struct();
for i = 1:n
  name_in = sprintf('subject%02d.mat',i);
  name_out = sprintf('s%02d_rep.mat',i);
  datapath_in = fullfile(datadir_in, name_in);
  datapath_out = fullfile(datadir_out, name_out);

  StagingContainer = load(datapath_in, data_varname_in);
  Data.raw = StagingContainer.(data_varname_in);
  clear StagingContainer;
  
  Replicated.visual = Data.raw(visual,:);
  Replicated.audio = Data.raw(audio,:);
  Replicated.semantic = Data.raw;
  clear Data;
  
  [~,filters(i).visual] = removeOutliers(Replicated.visual);
  [~,filters(i).audio] = removeOutliers(Replicated.audio);
  [~,filters(i).semantic] = removeOutliers(Replicated.semantic);
  
  save(datapath_out, '-struct', 'Replicated'); 
end

%% Compose metadata
metadata = struct();
for i = 1:n
  metadata(i).datacode = 'rep';
  metadata(i).subject = i;
  metadata(i).visual = visual;
  metadata(i).audio = audio;
  metadata(i).conceptcode = conceptcode;
  metadata(i).stimcode = stimcode;
  %
  metadata(i).filter(1).subset = 'visual';
  metadata(i).filter(1).label = 'rowfilter';
  metadata(i).filter(1).filter = filters(i).visual.words;
  metadata(i).filter(1).notes = 'Subset by modality, no aggregation.';
  %
  metadata(i).filter(2).subset = 'visual';
  metadata(i).filter(2).label = 'colfilter';
  metadata(i).filter(2).filter = filters(i).visual.voxels;
  metadata(i).filter(2).notes = 'Subset by modality, no aggregation.';
  %
  metadata(i).filter(3).subset = 'audio';
  metadata(i).filter(3).label = 'rowfilter';
  metadata(i).filter(3).filter = filters(i).audio.words;
  metadata(i).filter(3).notes = 'Subset by modality, no aggregation.';
  %
  metadata(i).filter(4).subset = 'audio';
  metadata(i).filter(4).label = 'colfilter';
  metadata(i).filter(4).filter = filters(i).audio.voxels;
  metadata(i).filter(4).notes = 'Subset by modality, no aggregation.';
  %
  metadata(i).filter(5).subset = 'semantic';
  metadata(i).filter(5).label = 'rowfilter';
  metadata(i).filter(5).filter = filters(i).semantic.words;
  metadata(i).filter(5).notes = 'No subsetting or aggregation.';
  %
  metadata(i).filter(6).subset = 'semantic';
  metadata(i).filter(6).label = 'colfilter';
  metadata(i).filter(6).filter = filters(i).semantic.voxels;
  metadata(i).filter(6).notes = 'No subsetting or aggregation.';

end

metapath_out = fullfile(datadir_out, 'metadata_rep.mat');
save(metapath_out, 'metadata');

%% Copy similarity structures into place
% They are already there, so I am not going to over engineer this...

%% Copy coordinates into place
% They are already there, so I am not going to over engineer this...

%% Pre-generate CV schemes
CV = defineCVBlocks(conceptcode, 'folds', 9, 'schemes', 11);
CV(:,12) = defineCVBlocks(conceptcode, 'LOO', true);