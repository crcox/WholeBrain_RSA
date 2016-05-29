# Simitar Toolbox
A release of the Botvinick lab at the Princeton Neuroscience Institute, written by francisco.pereira@gmail.com. This fork is maintained by Chris Cox at the University of Wisconsin-Madison, crcox@wisc.edu.

## Installation:

1. Add the toolbox directory to your MATLAB path ( `addpath(<directory>)` ). This will also work in Octave.

2. From MATLAB, go to the toolbox directory and try

```{matlab}
	mex findNeighbours.c
	mex simitar.c
	mex fastscoring.c
```

   If these work you will have much faster code for

   - preparing metas with createMetaFromMask
   - computing searchlight correlation matrices
   - computing similarity structure score maps

That's it! If you want to learn more please follow the [tutorial](http://www.princeton.edu/~fpereira/simitar/#tutorial) available on the [Simitar web page](http://www.princeton.edu/~fpereira/simitar).


## Data preparation:

In order to use any functions with your own data you will need to create a data
structure that maintains spatial information such as which voxels are neighbours
of which. There is a page on this topic linked to from the tutorial page, or
you can read [README.datapreparation.txt](README.datapreparation.txt), as it will explain how to create
this data structure from a binary mask indicating which voxels are in the brain.


## Functionality:

- Compute a searchlight similarity measure for the neighbourhood patch of each voxel
  function: [computeSimilarityMap.m](computeSimilarityMap.m)

- Scan for a particular similarity structure across all neighbourhood patches produced
  function: [computeSimilarityStructureMap.m](computeSimilarityStructureMap.m)

- Compute similarity measure between two sets of examples
  (euclidean, correlation, cosine and others)
  function: [computeExampleSimilarity.m](computeExampleSimilarity.m)

- Find the voxels that are in a given neighbourhood of each voxel
  (right now only supports cubic)
  function: [computeNeighboursWithinRadius.m](computeNeighboursWithinRadius.m)
