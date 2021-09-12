clear; clc; close all;

%% Adjacancy Matrix
A=[0 1 1 0 1 0 1 1 1 0 1 0
1 0 0 1 0 1 1 1 0 1 0 1
1 0 0 0 1 0 1 0 0 0 0 0
0 1 0 0 0 1 0 1 0 0 0 0
1 0 1 0 0 0 1 0 0 0 0 0
0 1 0 1 0 0 0 1 0 0 0 0
1 1 1 0 1 0 0 1 1 0 1 0
1 1 0 1 0 1 1 0 0 1 0 1
1 0 0 0 0 0 1 0 0 1 0 1
0 1 0 0 0 0 0 1 1 0 1 0
1 0 0 0 0 0 1 0 0 1 0 1
0 1 0 0 0 0 0 1 1 0 1 0];
%% vector of length n with the cluster index for each element
s=[1,1,3,3,3,3,1,1,2,2,2,2];
% Note that the P created is for a permuted version of A
[P,U,B] = clustered_SBD(A,s');





