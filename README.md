# Competing risk marginal frailty
This package is based on this paper AUTHORS

## Welcome
This package implements a semi-paramtric competing risk frailty model where each cluster has a vector of frailty parameters that are distributed multivariant normal.
Currently cluster sizes are constant but we accept missing cluster members.

## Model
The model accepts all observed failiure types and times and outputs the variance covariance matrix of the frailty vector, the estimated frailty for each cluster member and the cox regression that integrates the frailty.
The model expects clusters of the same size but its possible to indicate missing cluster members. The cox regression coefficient is different for every cluster member.
For example if several people take a test with three questions, then the coeffieicients will be different for each type question.

## example
The file example.pi contains a simple example of simulating data and fitting the model.

## performance
This is the simplest implementation. Parallel and GPU implementations exists and would be happily shared upon request.
