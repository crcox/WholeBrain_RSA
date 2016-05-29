% Computes similarity/distance measures within the searchlight for each voxel
% 
% As described in the tutorial, the code assumes two sets of examples -- A and B -- and their
% respective labels, and all the conditions present in one set should be present in the other.
% If you only have one set of examples, you should use it as both A and B.
%
% Prior to computing a similarity/distance matrix all examples of each condition in a set will
% be averaged together.
%
% The output, for each voxel, is a similarity/distance matrix between the average example of each
% condition in example set A and the average example of each condition in example set B.

% Input (mandatory):
%
% - measure - similarity/distance measure, can be: correlation | euclidean | cosine 
%
% - examples A - #examples A x #voxel
%   - each row can be an entire image or the voxels in a ROI
%   - each voxel can have raw signal or deconvolution coefficients
%
% - labels A   - #examples A x 1 (integers)
%   -numeric labels for the condition each example has

% - examples B - #examples B x #voxels
% - labels B - #examples B x 1 (integers)
%
% - neighbour information, in one of two forms
%   - 'neighbourInformation',<voxelsToNeighbours matrix>,<numberOfNeighbours matrix>
%       are such that:
%
%       <voxelsToNeighbours> - #voxels x maximum-#-neighbours
%       <numberOfNeighbours> - #voxels x 1 - actual-#-neighbours
%
%       to get the neighbours of voxel k (excluding itself)
%
%       neighbours = voxelsToNeighbours(k,1:numberOfNeighbours(k));
%
%   - 'meta',<meta structure> - the neighbour information is in the meta structure of a CMU IDM
%   format, overrides (or replaces) 'neighbourInformation' if available
%
%   - read README.datapreparation.txt or tutorial_data.html to see how to create this
%
% Output:
% - measurePerVoxel - #classes x #classes x #voxels
%   - measurePerVoxel(i,j,k) - measure of similarity/distance between the example of class <i>
%                              and the example of class <j> in voxel <k> searchlight
%   - the vectors of both classes contain voxel <k> and its adjacent neighbours
%
% Examples:
%
%   matrixPerLocation = computeSimilarityMap(measure,examples,labels,examples,labels,'meta',meta);
%
%
% History:
% 2013 Mar * - updated comments and documentation
% 2011 Sep * - cleaned up interface, added multiple sets of examples, different test types
% <wailing and gnashing of teeth>
% 2009 Feb 5 - fpereira@princeton.edu - created from earlier code
%
%   This file is part of Simitar
%
%   Simitar is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%   Simitar is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with Simitar.  If not, see <http://www.gnu.org/licenses/>.
%

function [measurePerVoxel] = computeSimilarityMap(varargin)

%
% process arguments and data
%
this = 'computeSimilarityMap';

if nargin > 0
  
  % first, figure out how many arguments before 'meta' or 'neighbourInformation'
  
  nargs = 0;
  idx = 1;
  while idx <= nargin
    argname = varargin{idx}; idx = idx + 1;
    if (isequal(argname,'meta') | isequal(argname,'neighbourInformation'))
      nargs = idx - 1; break;
    end
  end

  if nargs == 0
    fprintf('error: you need to specify a ''meta'' or a ''neighbourInformation'' parameter\n');return;
  end

  if nargs == 6
    nDatasets = 2;
    measure     = varargin{1};
    examples    = varargin{2};
    labels      = varargin{3};
    examples2   = varargin{4};
    labels2     = varargin{5};
  elseif nargs == 4
    nDatasets = 1;
    measure     = varargin{1};
    examples    = varargin{2};
    labels      = varargin{3};
  else
    fprintf('error: wrong number of mandatory arguments, either (measure,examples,labels,...) or (measure,examples1,labels1,examples2,labels2,...)\n');return;
  end
  
  switch measure
   case {'sloweuclidean'}
   case {'cosine'}
   case {'slowcorrelation'}
   case {'correlation','euclidean'}
    % requires the existence of a mex file, so test
    if exist('simitar.mexglx') | exist('simitar.mexa64') | exist('simitar.mexmac') | exist('simitar.mexmaci') | exist('simitar.mexmaci64') | exist('simitar.mexwin32') | exist('simitar.mexwin64')
      % OK
    else
      fprintf('error: using fastcorrelation/fasteuclidean requires compiling a C file with mex, try ''mex simitar.c''\n');return;
    end
    
   case {'L1'}
   otherwise
    fprintf('error: measure %s is not supported\n',measure);return;
  end
  
  [nExamples,nVoxels] = size(examples); n = nExamples; m = nVoxels;
  labelValues = unique(labels); nLabels = length(labelValues); k = nLabels;

  if nDatasets == 2
    % check that 
    [nExamples2,nVoxels2] = size(examples2);
    labelValues2 = unique(labels2);
    
    if (nVoxels2 ~= nVoxels) | ~isequal(labelValues2,labelValues)
      fprintf('error: the # of voxels and the label values must be the same across datasets\n');return;
    end
  end
  
  % defaults
  radius = 1;
  meta   = [];
  voxelsToNeighbours = [];
  numberOfNeighbours = [];
  labelsGroup = [];
  
  % process any additional arguments
  idx = nargs;
  while idx <= nargin
    argname = varargin{idx}; idx = idx + 1;
    switch argname
     case {'meta'}
      meta = varargin{idx}; idx = idx + 1;
     case {'neighbourInformation'}
      voxelsToNeighbours = varargin{idx}; idx = idx + 1;
      numberOfNeighbours = varargin{idx}; idx = idx + 1;
     otherwise
      whos
      pause
    
    end
  end

  % do some sanity checking
  
  if isempty(meta) & (isempty(voxelsToNeighbours) | isempty(numberOfNeighbours))
    fprintf('error: you must use either ''meta'' or ''neighbourInformation'' to specify neighbourhood information\n'); return;
  end
  
  if ~isempty(voxelsToNeighbours) & ~isempty(numberOfNeighbours)
    % trumps meta
  else
    % they must have been provided via meta
    if isempty(meta)
      % nope, do it for the entire image then, abort for now
      fprintf('error: %s requires neighbourhood information to be specified',this);return;
    else
      voxelsToNeighbours = meta.voxelsToNeighbours;
      numberOfNeighbours = meta.numberOfNeighbours;
    end
  end
  if ~isempty(voxelsToNeighbours)
    diameter = ceil(size(voxelsToNeighbours,2)^(1/3));
    radius   = (diameter-1)/2;
  end
else
  fprintf('syntax: see comments at the beginning of the function\n');return
end

%% figure out a few things

% which examples belong to which classes and groups
classes = unique(labels); nClasses = length(classes);

for c = 1:nClasses
  label = classes(c);
  maskClass{c} = (labels == label);
  indicesClass{c} = find(maskClass{c}); nPerClass(c) = length(indicesClass{c});
end
avgExampleClass = zeros(nClasses,nVoxels);

for c = 1:nClasses
  avgExampleClass(c,:) = mean(examples(indicesClass{c},:),1);
end

examples = avgExampleClass; n = nClasses;

if nDatasets == 2
  % we already know the # labels is the same (and they match)
  classes = unique(labels2); nClasses = length(classes);

  for c = 1:nClasses
    label = classes(c);
    maskClass{c} = (labels2 == label);
    indicesClass{c} = find(maskClass{c}); nPerClass(c) = length(indicesClass{c});
  end
  avgExampleClass = zeros(nClasses,nVoxels);

  for c = 1:nClasses
    avgExampleClass(c,:) = mean(examples2(indicesClass{c},:),1);
  end

  examples2 = avgExampleClass;
end

%
% Compute distances
%

fprintf('%s: computing distances\n',this);

measurePerVoxel = zeros(nClasses,nClasses,nVoxels);

switch measure
  
 case {'correlation'}
  % code will output into measurePerVoxel, to avoid having to reallocate large amounts of memory
  if nDatasets == 2
    simitar(examples',examples2',radius,voxelsToNeighbours',numberOfNeighbours',measurePerVoxel,1);
  else
    simitar(examples',examples',radius,voxelsToNeighbours',numberOfNeighbours',measurePerVoxel,1);
  end
    
 case {'euclidean'}
  % code will output into measurePerVoxel, to avoid having to reallocate large amounts of memory
  if nDatasets == 2
    simitar(examples',examples2',radius,voxelsToNeighbours',numberOfNeighbours',measurePerVoxel,2);
  else
    simitar(examples',examples',radius,voxelsToNeighbours',numberOfNeighbours',measurePerVoxel,2);
  end
  
 otherwise
  
  if nDatasets == 1

    % 1 set of examples
    len = nClasses;

    for v = 1:nVoxels
      neighbours = [v,voxelsToNeighbours(v,1:numberOfNeighbours(v))];

      % prepare a normalized version of the dataset
      switch measure
       case {'slowcorrelation'}
        tmp = (examples(:,neighbours) - repmat(mean( examples(:,neighbours),2),1,numberOfNeighbours(v)+1))./ repmat(std( examples(:,neighbours),0,2),1,numberOfNeighbours(v)+1);      
       otherwise
      end
      
      for c = 1:nClasses
        switch measure
         case {'sloweuclidean'}
          measurePerVoxel(:,c,v) = ...
              sum((repmat(examples(c,neighbours),len,1)-examples(:,neighbours)).^2,2);
         case {'L1'}
          measurePerVoxel(:,c,v) = ...
              sum(abs(repmat(examples(c,neighbours),len,1)-examples(:,neighbours)),2);
          
         case {'slowcorrelation'}
          measurePerVoxel(:,c,v) = sum(repmat(tmp(c,:),len,1) .* tmp,2)/numberOfNeighbours(v);
          
          
         case {'cosine'}
          measurePerVoxel(:,c,v) = sum(...
              repmat(examples(c,neighbours),len,1) .* examples(:,neighbours), 2) ./ ...
              sqrt(sum((examples(c,neighbours)).^2,2)*sum((examples(:,neighbours)).^2,2));
        end
      end
      
      if ~rem(v,1000); fprintf('\t%d',v); end
    end

  else
    
    % 2 sets of examples
    
    len = nClasses;
    
    for v = 1:nVoxels
      %  for v = 91800:nVoxels
      neighbours = [v,voxelsToNeighbours(v,1:numberOfNeighbours(v))];
      
      % prepare a normalized version of the dataset
      switch measure
       case {'slowcorrelation'}
        tmp1 = (examples(:,neighbours) - repmat(mean( examples(:,neighbours),2),1,numberOfNeighbours(v)+1))./ repmat(std( examples(:,neighbours),0,2),1,numberOfNeighbours(v)+1);      
        tmp2 = (examples2(:,neighbours) - repmat(mean( examples2(:,neighbours),2),1,numberOfNeighbours(v)+1))./ repmat(std( examples2(:,neighbours),0,2),1,numberOfNeighbours(v)+1);      
       otherwise
      end
      
      for c = 1:nClasses
        switch measure
         case {'sloweuclidean'}
          measurePerVoxel(:,c,v) = ...
              sum((repmat(examples(c,neighbours),len,1)-examples(:,neighbours)).^2,2);
         case {'L1'}
          measurePerVoxel(:,c,v) = ...
              sum(abs(repmat(examples(c,neighbours),len,1)-examples(:,neighbours)),2);
         case {'slowcorrelation'}
          measurePerVoxel(c,:,v) = sum(repmat(tmp1(c,:),len,1) .* tmp2,2)/numberOfNeighbours(v);
          
         case {'cosine'}
          measurePerVoxel(:,c,v) = sum(...
              repmat(examples(c,neighbours),len,1) .* examples(:,neighbours), 2) ./ ...
              sqrt(sum((examples(c,neighbours)).^2,2)*sum((examples(:,neighbours)).^2,2));
        end
      end

      %    fprintf('%d %d %d %d %d\n',v,size(measurePerVoxel),nClasses);
      
      %clf;
      %imagesc(measurePerVoxel(:,:,v)); axis square;
      %pause
      
      if ~rem(v,1000); fprintf('\t%d',v); end
    end
    
    
  end

end
  
fprintf('\n');
