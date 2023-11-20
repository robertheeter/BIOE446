close all
clc

%% Problem 1: 1D pattern formation in the Gierer-Meinhardt model
disp('PROBLEM 1')
clear

%% Problem 1, Part A
P = [0.1 1.25 0.07];

tmax = 100;
n = 100;
L = 10;

soln = GM_solve(P , tmax, n, L, [1 2 3]);

%% Problem 1, Part B
rxn_dfsn_gs(0.09, 0.059, 1, 0.5, InitState='uSpots');

%% Problem 1, Part C


%% Problem 2: 2D pattern formation in the Gray-Scott model
disp('PROBLEM 2')
clear

%% Problem 2, Part A


%% Problem 2, Part B

%% Problem 2, Part C


%% Functions

% None
