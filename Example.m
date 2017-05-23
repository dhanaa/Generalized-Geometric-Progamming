%% example file for testing Generalized Geometric Programming problem
% Duo, Hong Kong University of Science and Technology
% dhanaa@connect.ust.hk 
% 2017/05/22

%% write the original problem into matrix form
% use vector to represent monomial
% use struct to represent each line, one for positive, one for
% negative
% example in page 20: % http://maranas.che.psu.edu/pub/1997/
%Maranas_and_Floudas,_Computers_and_Chem._Eng.,_1997.pdf

% 2 variables, 
% min x
% st  0.25 * x + 0.5 * y - (1/16) * x^2 - (1/16)*y^2 -1 <= 0
%     (1/14) * x^2 + (1/14) * y^2 + 1 -(3/7) * x - (3/7) * y <=0
%      1 <= x <= 5.5
%      1 <= y <= 5.5

close all; clear all; clc;

% objective function
obj.positive = [1 1 0];
obj.negative = 0;

% constraint
con1.positive = [0.25 1 0; 0.5 0 1];
con1.negative = [1/16 2 0; 1/16 0 2; 1 0 0];

con2.positive = [1/14 2 0; 1/14 0 2; 1 0 0];
con2.negative = [3/7 1 0; 3/7 0 1];

constraint = [con1, con2]; % stack all constraints

lowerbound = [1, 1]; % global bounds for variables
upperbound = [5.5, 5.5];

N = 2; % number of variables

[obj, var] = GGPSolver(N, obj,constraint,lowerbound,upperbound);
