function [results, params, varargout] = LoadAndProcessResults(resultsdir, csvfilename, mridir, opts, varargin)
  % Check function arguments
  p = inputParser();
  addRequired(p,'resultsdir',@ischar);
  addRequired(p,'csvfilename',@ischar);
  addRequired(p,'mridir',@ischar);
  addRequired(p,'opts',@isstruct);
  addParameter(p, 'PermTest', false, @islogical);
  parse(p,resultsdir, csvfilename, mridir, opts, varargin{:});

  resultsdir = p.Results.resultsdir;
  csvfilename = p.Results.csvfilename;
  mridir = p.Results.mridir;
  PermTest = p.Results.PermTest;
  fields = fieldnames(p.Results.opts);
  values = struct2cell(p.Results.opts);
  optscell = [fields(:)';values(:)'];

  % Check opts
  q = inputParser();
  addParameter(q,'resultsfile',@ischar)
  addParameter(q,'paramsfile',@ischar)
  addParameter(q,'datadir',@ischar)
  addParameter(q,'datafilefmt',@ischar)
  addParameter(q,'datavarname',@ischar)
  addParameter(q,'metafile',@ischar)
  addParameter(q,'metavarname',@ischar)
  addParameter(q,'coordsfile',@ischar)
  addParameter(q,'coordsvarname',@ischar)
  addParameter(q,'coordslabel',@ischar)
  addParameter(q,'writefields',@iscellstr)
  parse(q,optscell{:});

  opts = q.Results;

  % Load metadata
  tmp = load(fullfile(opts.datadir, opts.metafile), opts.metavarname);
  metadata = tmp.(opts.metavarname);
  tmp = load(fullfile(opts.datadir, opts.coordsfile), opts.coordsvarname);
  coords = tmp.(opts.coordsvarname);

  % Generate filenames
  [p,f,e] = fileparts(csvfilename);
  matfilename = fullfile(p,strcat(f,'.mat'));
  permfilename = fullfile(p,strcat(f,'_aggregated.mat'));

  if exist(matfilename, 'file')
    load(matfilename, 'results','params');
    if ~exist(csvfilename, 'file')
      % Write Results to csv (if needed)
      WriteTable(csvfilename, results, params, 'fields', opts.writefields, 'overwrite', false);
    end
  else
    % Load Results
    [results, params] = LoadResults( ...
      'ResultDir', resultsdir, ...
      'DataDir', opts.datadir, ...
      'MetadataFile', opts.metafile, ...
      'MetadataVarname', opts.metavarname, ...
      'ResultFile', opts.resultsfile, ...
      'ParamFile', opts.paramsfile);

    % Save results
    save(matfilename, 'results','params');

    % Write Results to csv
    WriteTable(csvfilename, results, params, 'fields', opts.writefields, 'overwrite', false);
  end

  % Aggregate over permutations
  if PermTest
    if exist(permfilename, 'file')
      load(permfilename, 'results_perm');
    else
      results_perm = AggregatePermutationResults(results);
      save(permfilename, 'results_perm');
    end

    % Pair coordinates with discovered voxels
    results_perm = PairWithCoordinates(results_perm, metadata, coords, opts.coordslabel);
    varargout{1} = results_perm;

    warning('off','MATLAB:MKDIR:DirectoryExists');
    mkdir(mridir)
    mkdir(fullfile(mridir,'count'))
    mkdir(fullfile(mridir,'count','txt'));
    mkdir(fullfile(mridir,'count','csv'));
    mkdir(fullfile(mridir,'count','figures'));
    mkdir(fullfile(mridir,'nodestrength'))
    mkdir(fullfile(mridir,'nodestrength','txt'));
    mkdir(fullfile(mridir,'nodestrength','csv'));
    mkdir(fullfile(mridir,'nodestrength','figures'));
    warning('on','MATLAB:MKDIR:DirectoryExists');
    WriteVoxels(fullfile(mridir,'count','txt'), results_perm, 'CoordsLabel', opts.coordslabel, 'VoxelValue', 'count_nz_rows');
    WriteVoxels(fullfile(mridir,'nodestrength','txt'), results_perm, 'CoordsLabel', opts.coordslabel, 'VoxelValue', 'mean_nodestrength');
  else
    % Pair coordinates with discovered voxels
    results = PairWithCoordinates(results, metadata, coords, opts.coordslabel);

    % Write voxel data to text
    warning('off','MATLAB:MKDIR:DirectoryExists');
    mkdir(mridir)
    mkdir(fullfile(mridir,'txt'));
    mkdir(fullfile(mridir,'csv'));
    mkdir(fullfile(mridir,'figures'));
    warning('on','MATLAB:MKDIR:DirectoryExists');
    WriteVoxels(fullfile(mridir,'txt'), results, 'CoordsLabel', opts.coordslabel);
end
