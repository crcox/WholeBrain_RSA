% Load and process final/plotting data
%
%% Add paths
addpath('~/src/WholeBrain_RSA/util');
addpath('~/src/WholeBrain_RSA/dependencies/jsonlab/');

%% Load data
rdir  = '~/MRI/Manchester/results/WholeBrain_RSA/semantic/audio/grOWL2/final/plotting/';
ddir  = '~/MRI/Manchester/data/avg';
mfile = 'metadata_avg.mat';
mvar  = 'metadata';
load(fullfile(ddir,mfile), mvar);

rfile = 'results.mat';
pfile = 'params.json';
skip = {};
[results, params] = LoadResults('ResultDir', rdir, 'DataDir', ddir, ...
    'MetadataFile', mfile, 'MetadataVarname', mvar, ...
    'ResultFile', rfile, 'ParamFile', pfile, ...
    'SkipFields', skip);

%% Inspect
results(1)

%% Add coordinates to the results object
cfile = fullfile(ddir,'coords_avg_adj.mat');
cvar = 'coords';
load(cfile, cvar);
for i = 1:numel(results)
    R = results(i);
    s = R.subject;
    z = strcmp('colfilter', {metadata(s).filter.label});
    v = metadata(s).filter(z).filter;
    w = R.nz_rows;
    if R.bias
        w = w(1:(end-1));
    end
    xyz = coords(s).mni(v,:);
    xyz = xyz(w,:);
    label = 'mni';
    results(i).coords(1).label = label;
    results(i).coords(1).xyz = xyz;
end

%% Add node strength to the results object
for i = 1:numel(results)
    R = results(i);
    s = R.subject;
    w = R.nz_rows;
    Uz = R.Uz;
    if R.bias
        w = w(1:(end-1));
        Uz = Uz(1:(end-1),:);
    end
    Uz = Uz(w,:);
    Wz = Uz * Uz';
    results(i).nodestrength = sum(abs(Wz))/max(sum(abs(Wz)));
end

%% Dump
tdir  = fullfile(rdir,'solutionmaps','txt');
mkdir(tdir);
WriteVoxels(tdir,results,'VoxelValue','nodestrength')

%% Package for Tim
mkdir(fullfile(rdir,'Tim'));
for i = 1:numel(results)
    R = results(i);
    s = R.subject;
    w = R.nz_rows;
    Uz = R.Uz;
    if R.bias
        w = w(1:(end-1));
        Uz = Uz(1:(end-1),:);
    end
    Uz = Uz(w,:);
    xyz = R.coords(1).xyz;
    uname = sprintf('%02d_Uz.csv',i);
    cname = sprintf('%02d_mni.csv',i);
    upath = fullfile(rdir,'Tim',uname);
    cpath = fullfile(rdir,'Tim',cname);
    csvwrite(upath, Uz);
    csvwrite(cpath, xyz);
end