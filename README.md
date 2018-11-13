# MasterThesis

This repository algorithms I wrote during my Master Thesis at Department of Physics of University of Warsaw. 
It contains implementation of (1) for cavity with displaced mirror: except for constant parameters, input is just initial state of light which goes through a quantum channel. Then, it's time evolution and quantum Fisher information (QFI) is calculated.

Independent implementations in Mathematica and Matlab (cvxr (2), because iterative algorithm for optimization solves a convex problem) contain analysis of state in time: its coefficients, QFI and dependency on some parameters like coupling constant and squeezing angle.

TODO:
* add richer description to readme
* add link to current version of notes in Overleaf
* fix Matlab code so that it works well for all t, g an dimension, not only for some ranges.

(1) K. Macieszczak. Quantum Fisher Information: Variational principle and simpleiterative algorithm for its efficient computation. arXiv:1312.1356v1, 2013
(2) CVX:  Matlab  Software  for  Disciplined  Convex  Programming. http://cvxr.com/cvx/24



