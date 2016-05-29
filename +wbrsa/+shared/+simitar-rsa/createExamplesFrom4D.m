%
% Transforms a dataset in NIfTI format into examples+meta
%
% input:
% - data - 4D array
% - meta - meta structure obtained with createMetaFromMask
%
% output: (#voxels is either whole volume or those in brain mask, if provided)
% - examples - #volumes x #voxels
%
% Notes:
% - depends on having the FSL MATLAB functions in the path
%

function [examples] = createExamplesFrom4D(varargin)

%
% process parameter
% 

if nargin < 2
  fprintf('syntax: createExamplesAndMetaFrom4D(<dataset in 4D>,<meta>\n');
  return;
end

data = varargin{1}; dims = size(data);
meta = varargin{2}; mdims = meta.dimensions;

if ~isequal(mdims(1:3),dims(1:3))
  fprintf('error: dataset and mask 3D dimensions are different!\n');return;
end

%
% create examples from data and meta
%

m = length(meta.indicesIn3D);
n = size(data,4);

examples = zeros(n,m);
volume   = zeros(meta.dimensions);

for t = 1:n
  volume = data(:,:,:,t);
  examples(t,:) = volume(meta.indicesIn3D);
end
