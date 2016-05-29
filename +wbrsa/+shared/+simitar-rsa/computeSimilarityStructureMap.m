% Computes similarity measure within the searchlight for each voxel
% and scores it for a given structure, using permutation tests to obtain a p-value
%
% (for more details about the examples and labels arguments please refer to the help
%  for the computeSimilarityMap.m function)
%
% Input (mandatory):
%
% - measure - similarity/distance measure, can be: correlation | euclidean | cosine (slow)
%
% - examples A - #examples A x #voxel
% - labels A   - #examples A x 1 (integers)
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
%   read README.datapreparation.txt or tutorial_data.html to see how to create this
%
% - 'structureMatrix',<ternary #conditions x #conditions>
%
%   A ternary matrix specifying how the values in each similarity/distance matrix contribute
%   to a numeric score that should be high for matrices that have the desired structure
%   (please see examples below or read the tutorial for more details).
%
%   The matrix gets multiplied elementwise by each similarity/distance matrix, and all the
%   entries of the result are summed to produce the numeric score.
%
% Input (optional):
%
% - statistical test information, in one of these forms
%   - 'permutationTest',<test type>,<# permutations>
%
%      where test type is either 'overExamples' (preferable) or 'overMatrices'
%      (please see tutorial for more details)
%
% - example group information, as
%   - 'groupLabels',<labels group A>,<labels group B>
%
%   - this is used for permutation tests, as label permutations will happen within group only
%   - a run is the typical grouping unit
%
% Output:
%
% - structureMap - 1 x #voxels
%   - each voxel will contain the score obtained by multiplying the structure score matrix by
%     the respective similarity/distance matrix, elementwise, and summing all entries of the result
%
% - pvalueMap    - 1 x #voxels
%   - if permutation tests were requested, the p-value of the structure score obtained at that
%   voxel under the permutation distribution specified (or [] if no tests were requested).
%
%
% Examples:
%
% using this example structure matrix, where
% condition 1 is similar to condition 4 but dissimilar to conditions 2 and 3
%
%  structureMatrix = [[0 -1 -1 1];[-1 0 0 -1];[-1 0 0 -1];[1 -1 -1 0]];
%
% this will produce a map where each voxel has a value indicating whether this structure is present
%
% [structureScoreMap] = computeSimilarityStructureMap('correlation',examples,labels,examples,labels,'meta',meta,'similarityStructure',structureMatrix);
%
% if you want a permutation test with 1000 permutations over all examples 
% 
% testType = 'overExamples';
% nPermutations = 1000;
%      [structureScoreMap,structurePvalueMap] = computeSimilarityStructureMap('correlation',examples,labels,examples,labels,'meta',meta,'similarityStructure',structureMatrix,'permutationTest',testType,nPermutations,'groupLabels',labelsRun,labelsRun);
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

function [structureScoreMap,structurePvalueMap] = computeSimilarityStructureMap(varargin)

%
% process arguments and data
%
this = 'computeSimilarityStructureMap';
structureScoreMap = []; structurePvalueMap = [];

if nargin < 9
  fprintf('error: wrong number of mandatory arguments, (measure,examples1,labels1,examples2,labels2,''meta'',<meta>,''similarityStructure'',<structure matrix>,...)\n');return;  
else
  %% process first batch of mandatory arguments
  
  measure     = varargin{1};
  examplesA   = varargin{2};
  labelsA     = varargin{3};
  examplesB   = varargin{4};
  labelsB     = varargin{5};
  
  % check them
  
  switch measure
   case {'sloweuclidean'}
   case {'cosine'}
   case {'slowcorrelation'}
   case {'correlation','euclidean'}
    % requires the existence of a mex file, so test
    if exist('simitar.mexglx') | exist('simitar.mexa64') | exist('simitar.mexmac') | exist('simitar.mexmaci') ...
          | exist('simitar.mexmaci64') | exist('simitar.mexwin32') | exist('simitar.mexwin64')
      % OK
    else
      fprintf('error: using fastcorrelation/fasteuclidean requires compiling a C file with mex, try ''mex simitar.c''\n');return;
    end
    
   case {'L1'}
   otherwise
    fprintf('error: measure %s is not supported\n',measure);return;
  end
  
  % get info about both sets of examples and check that they match
  
  [nExamplesA,nVoxelsA] = size(examplesA); n = nExamplesA; m = nVoxelsA;
  labelValuesA = unique(labelsA); nLabelsA = length(labelValuesA); 
  
  if isequal(examplesA,examplesB)
    datasetsAreTheSame = 1;
    nExamplesB = nExamplesA;
    nVoxelsB = nVoxelsA;
    labelValuesB = labelValuesA;
  
  else
    datasetsAreTheSame = 0;
  
    [nExamplesB,nVoxelsB] = size(examplesB);
    labelValuesB = unique(labelsB);
    
    if (nVoxelsB ~= nVoxelsA) | ~isequal(labelValuesB,labelValuesA)
      fprintf('error: the # of voxels and the label values must be the same across datasets\n');return;
    end
  end

  % now that we know they match
  nVoxels = nVoxelsA;
  labelValues = labelValuesA; classes = labelValues;
  nLabels = nLabelsA; nClasses = nLabels;
  k = nLabels;
  m = nVoxels;
    
  %% process additional arguments

  radius = 1;
  meta   = [];
  voxelsToNeighbours = [];
  numberOfNeighbours = [];
  labelsGroup = [];
  nperms = 0;
  labelsGroupA = [];
  labelsGroupB = [];
  
  % process any additional arguments
  idx = 6;
  while idx <= nargin
    argname = varargin{idx}; idx = idx + 1;
    switch argname
     case {'similarityStructure'}
      similarityStructure = varargin{idx}; idx = idx + 1;
     case {'meta'}
      meta = varargin{idx}; idx = idx + 1;
     case {'neighbourInformation'}
      voxelsToNeighbours = varargin{idx}; idx = idx + 1;
      numberOfNeighbours = varargin{idx}; idx = idx + 1;
     case {'permutationTest'}     
      permutationTestType = varargin{idx}; idx = idx + 1;
      nperms = varargin{idx}; idx = idx + 1;
      switch permutationTestType
       case {'overMatrices','overExamples'}
       otherwise
        fprintf('error: arg should be ''permutationTest'',<test type>,<#perms>\n');return;
      end

     case {'groupLabels'}
      labelsGroupA = varargin{idx}; idx = idx + 1;
      labelsGroupB = varargin{idx}; idx = idx + 1;
     otherwise
      fprintf('error: unknown argument %s\n',argname);return;
    end
  end
  
  % and check
  
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

  if isempty(similarityStructure)
    fprintf('error: you must specify a similarityStructure'); return;
  end
  
  [k1,k2] = size(similarityStructure);
  if (k~=k1) | (k~=k2)
    fprintf('error: the similarity structure matrix must have the same dimensions as the matrix for each voxel.\n');return;
  end
  
  if ~isequal(similarityStructure,similarityStructure')
    fprintf('error: the similarity structure matrix must be symmetric.\n'); return;
  end
  
  % check the statistics are valid ones and do statistic-specific preprocessing
  
  indicesDiag = sub2ind([k k],(1:k)',(1:k)');
  indicesOne   = find(similarityStructure(:)== 1);  nOne  = length(indicesOne);
  indicesZero  = find(similarityStructure(:)== 0);  nZero = length(indicesZero);
  indicesMinus = find(similarityStructure(:)==-1); nMinus = length(indicesMinus);

  if (nZero+nOne+nMinus) ~= (k*k)
    fprintf('warning: your similarity structure score matrix is not ternary (-1,0,1).\n');
    useTernaryMatrix = 0;
    %    fprintf('error: the similarity structure matrix must be ternary (-1,0,1).\n'); return;
  else
    fprintf('using a ternary similarity structure score matrix\n');
    useTernaryMatrix = 1;
  end
 
  
  % figure out whether the MEX version exists

  if exist('fastscoring.mexglx') | exist('fastscoring.mexa64') | exist('fastscoring.mexmac') | exist('fastscoring.mexmaci') | exist('fastscoring.mexmaci64') | exist('fastscoring.mexwin32') | exist('fastscoring.mexwin64')
    useMEX = 1;
  else
    useMEX = 0;
  end
end

%% figure out a few things

%% which examples belong to which classes and groups, average example per class

for c = 1:nClasses
  label = classes(c);
  maskClassA{c} = (labelsA == label);
  indicesClassA{c} = find(maskClassA{c}); nPerClassA(c) = length(indicesClassA{c});
end

% check whether the permutation test type selected is compatible with the #examples
if nperms > 0
  switch permutationTestType
   case {'overExamples'}
    if sum(nPerClassA<repmat(5,1,nClasses))
      fprintf('error: too few examples to do permutation test over examples, please select a test over matrices\n');return;
    end
   case {'overMatrices'}
    % test s are done inside the branch for this
  end
end


averageExampleClassA = zeros(nClasses,nVoxels);
for c = 1:nClasses
  averageExampleClassA(c,:) = mean(examplesA(indicesClassA{c},:),1);
end

if datasetsAreTheSame
  maskClassB = maskClassA;
  indicesClassB = indicesClassA;
  nPerClassB = nPerClassA;
  averageExampleClassB = averageExampleClassA;
else
  for c = 1:nClasses
    label = classes(c);
    maskClassB{c} = (labelsB == label);
    indicesClassB{c} = find(maskClassB{c}); nPerClassB(c) = length(indicesClassB{c});
  end

  averageExampleClassB = zeros(nClasses,nVoxels);
  for c = 1:nClasses
    averageExampleClassB(c,:) = mean(examplesB(indicesClassB{c},:),1);
  end
end

%examples = averageExampleClass; n = nClasses;
%examples2 = averageExampleClass;

%
% Compute distances (a kind of cross-validation)
%

fprintf('%s: computing distances\n',this);

% will be used for both the main computation and the permutation test
measurePerVoxel = zeros(nClasses,nClasses,nVoxels);

switch measure
 case {'correlation'}
  % code will output into measurePerVoxel, to avoid having to reallocate large amounts of memory
  simitar(averageExampleClassA',averageExampleClassB',radius,voxelsToNeighbours',numberOfNeighbours',measurePerVoxel,1);

 case {'euclidean'}
  % code will output into measurePerVoxel, to avoid having to reallocate large amounts of memory
  simitar(averageExampleClassA',averageExampleClassB',radius,voxelsToNeighbours',numberOfNeighbours',measurePerVoxel,2);
  
 otherwise
  examples  = averageExampleClassA;
  examples2 = averageExampleClassB;
  
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
    %if ~rem(v,1000); fprintf('\t%d',v); end
  end
end
  
%
% compute score
% 

structureScoreMap  = zeros(1,m);
structurePvalueMap = ones(1,m);

% used for fast calculation, balance the weight of +1 and -1
tmps = similarityStructure;

if useTernaryMatrix
  tmps(indicesOne)   = tmps(indicesOne)/nOne;
  tmps(indicesMinus) = tmps(indicesMinus)/nMinus;
else
  % user takes care of balancing
end
  
  
tic;
if useMEX
  fastscoring(measurePerVoxel(:),tmps(:),k,m,structureScoreMap);
else
  tmp = reshape(measurePerVoxel,(k*k),m);

  for v = 1:m
    structureScoreMap(v) = sum(tmp(:,v) .* tmps(:));
  end
end
t = toc;

%
% permutation tests
%

if nperms > 0

  fprintf('computing p-values using a %d-permutation test over ',nperms);
  
  switch permutationTestType
    
   case {'overExamples'}
    fprintf('examples\n');
    fprintf('#permutations: ');
    
    nperms = nperms - 1; % leave 1 for the true labels
    structureScoreMapPermuted  = zeros(1,m);
    structurePvalueMap = zeros(1,m);structurePvalueMap = zeros(1,m);
    
    %% figure out a few things
    
    if isempty(labelsGroupA) | isempty(labelsGroupB)
      % just assume all the examples are in one large group
      labelsGroupA = ones(nExamplesA,1);
      labelsGroupB = ones(nExamplesB,1);
    end
    
    groupsA = unique(labelsGroupA); nGroupsA = length(groupsA);
    for ig = 1:nGroupsA
      indicesGroupA{ig} = find(labelsGroupA == groupsA(ig));
    end

    groupsB = unique(labelsGroupB); nGroupsB = length(groupsB);
    for ig = 1:nGroupsB
      indicesGroupB{ig} = find(labelsGroupB == groupsB(ig));
    end
    
    %% main loop
    
    indicesA = (1:nExamplesA)';
    indicesB = (1:nExamplesB)';
    
    if nperms < 1000; interval = 100; else; interval = 100; end

    tstart = tic;
    
    for ip = 1:nperms
      if ~rem(ip,interval); fprintf(' %d',ip); end

      if ip == 10
        telapsed = toc(tstart);
        nloops = ceil(nperms/10 -1);
        fprintf('estimated time to completion: %d second(s)\n',ceil(nloops*telapsed));
      end
        
      % permute indices inside each group
      
      indicesApermuted = zeros(nExamplesA,1);
      indicesBpermuted = zeros(nExamplesB,1);
      
      for ig = 1:nGroupsA
        indicesApermuted(indicesGroupA{ig}) = indicesA(indicesGroupA{ig}(randperm(length(indicesGroupA{ig}))));
      end
      
      if datasetsAreTheSame
        indicesBpermuted = indicesApermuted;
      else
        for ig = 1:nGroupsB
          indicesBpermuted(indicesGroupB{ig}) = indicesB(indicesGroupB{ig}(randperm(length(indicesGroupB{ig}))));
        end
      end

      %disp([labelsGroupA,labelsGroupB,indicesA,indicesB]);pause
      %disp([labelsGroupA,indicesA,indicesApermuted,indicesBpermuted]);pause
      
      % compute the average example for each class with labels permuted

      for c = 1:nClasses
        label = classes(c);
        indicesClassA{c} = find(labelsA(indicesApermuted) == label);
      end
      averageExampleClassA = zeros(nClasses,nVoxels);
      for c = 1:nClasses
        averageExampleClassA(c,:) = mean(examplesA(indicesClassA{c},:),1);
      end

      for c = 1:nClasses
        label = classes(c);
        indicesClassB{c} = find(labelsB(indicesBpermuted) == label);
      end
      averageExampleClassB = zeros(nClasses,nVoxels);
      for c = 1:nClasses
        averageExampleClassB(c,:) = mean(examplesB(indicesClassB{c},:),1);
      end

      % compute similarity and score
      
      switch measure
       case {'correlation'}
        % code will output into measurePerVoxel, to avoid having to reallocate large amounts of memory
        simitar(averageExampleClassA',averageExampleClassB',radius,voxelsToNeighbours',numberOfNeighbours',measurePerVoxel,1);
        
       case {'euclidean'}
        % code will output into measurePerVoxel, to avoid having to reallocate large amounts of memory
        simitar(averageExampleClassA',averageExampleClassB',radius,voxelsToNeighbours',numberOfNeighbours',measurePerVoxel,2);
        
       otherwise
        examples  = averageExampleClassA;
        examples2 = averageExampleClassB;        
        len = nClasses;
        
        for v = 1:nVoxels
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
          %if ~rem(v,1000); fprintf('\t%d',v); end
        end
      end; % switch on measure

      if useMEX
        fastscoring(measurePerVoxel(:),tmps(:),k,m,structureScoreMapPermuted);
      else
        tmp = reshape(measurePerVoxel,(k*k),m);
        
        for v = 1:m
          structureScoreMapPermuted(v) = sum(tmp(:,v) .* tmps(:));
        end
      end
      
      % count
      
      switch measure
       case {'euclidean','sloweuclidean','correlation','slowcorrelation','cosine'}
        % measures where higher is better
        structurePvalueMap = structurePvalueMap + (structureScoreMapPermuted>=structureScoreMap);
        % no longer used (just use inverse sign matrices)
       %case {'sloweuclidean','L1'}
        % measures where lower is better
        %structurePvalueMap = structurePvalueMap + (structureScoreMapPermuted<=structureScoreMap);
      end
      
    end; % for over #permutations
    fprintf('\n');
    
    structurePvalueMap = (structurePvalueMap+1)/ nperms;
    
    
   case {'overMatrices'}
    fprintf('matrices\n');
    
    if nperms
      counts = zeros(1,m);
      npairs = k*(k-1)/2;

      if k <= 3
        fprintf('warning: not enough condition pairs to run a matrix permutation test\n');
      else
        tmpsv = transformMatrixToLT(tmps,0);
        results = zeros(1,m);
        
        if k == 4
          fprintf('warning: 4 conditions, will do all possible permutations (720)\n');
          permutations = perms(1:6);
          nperms = 720;
        else
          fprintf('computing permutation test p-values with %d permutations\n',nperms);
        end

        %fprintf('Estimated time %1.2f seconds\t',nperms*t);

        for ip = 1:nperms
          if k == 4
            rp = permutations(ip,:);
          else
            rp = randperm(npairs);
          end
          rtmp = tmpsv(rp);
          stmp = transformLTtoMatrix(rtmp,k,0);
          stmp = stmp + stmp';
          
          %  results(v) = sum(tmp(:,v) .* stmp(:));
          fastscoring(measurePerVoxel(:),stmp(:),k,m,results);
          
          counts = counts + (results>=structureScoreMap);
        end

        %t=toc;
        %fprintf('took %1.2f seconds\n',t);

        if k == 4
          counts = counts / nperms;
        else
          % add 1 to count the true labels as having occurred once
          % (conservative if they have already occurred, but faster not to check every time)
          counts = (1+counts) / nperms;
        end
        
        %    imagesc([counts;structureScoreMap]);pause
        %    clf;
        %similarityStructure
        %dimx = sqrt(m);
        %clf;
        %subplot(1,2,1); imagesc(reshape(1-counts,[dimx dimx]),[0 1]); colorbar('vert');axis square;
        %subplot(1,2,2); imagesc(reshape(structureScoreMap,[dimx dimx]),[0 1]); colorbar('vert');axis square;
        %pause
      end
      
      structurePvalueMap = counts;
    end; % if nperms
    
  end; % of switch on permutation type
  
  
end; % if nperms > 0
