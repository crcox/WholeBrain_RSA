# Whole Brain RSA
The objective of this procedure is to learn how a system encodes a particular similarity structure. For example, let's sat you have a symmetric matrix `S` that expresses the similarity among `n` items---for example, how similar each of a set of words are in their meaning---and you want to learn how some system---for example, the brain---represents that similarity structure. Let the brain data be represented by an `n by d` matrix `X`, where `n` is the number of items and `d` is the number of voxels it takes to represent the brain activity for each item at some spatial resolution. The goal is to find a sparse subset of the `d` voxels that encode the similarity structure in `S`. That is, we want to find the best sparse set of voxels, `v*`, so that `corr(S, X_{v*} * U * X_{v*}^T)` is as high as possible, where `U` is a matrix that weights the relationships among voxels and each voxels overall importance. Conceptually, this procedure attempts to find `v*` and solve for `U`.

# Instructions
Because this is a computationally expensive operation, with a free parameter to be tuned, this code has been written with an eye towards distributed computing. That is, the main function `WholeBrain_RSA` is written so that it expects to find a file called `params.yaml` in the current working directory. This file should contain all of the parameters and specifications required for loading loading and processing data, and for generally instructing the algorithm on how to perform. 

The process is abstracted into three main functions:
## WholeBrain_RSA
`WholeBrain_RSA` is the outter-most wrapper, where data are loaded and processed. This is where the `yaml` file will be parsed into Matlab.

## learn_similarty_encoding
`learn_similarity_encoding` is responsible for the more technical transformations of the target similarity structure and the data.

## Adlas1
`Adlas1` actually executes the group_lasso algorithm.

