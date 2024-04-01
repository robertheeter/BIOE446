close all
clc

%% Problem 1: 1D pattern formation in the Gierer-Meinhardt model
disp('PROBLEM 1')
clear

%% Problem 1, Part A
close all
clc
clear

tmax = 100;
n = 100;
L = 10;

% P1 = [0.1 1.25 0.07];
% GM_solve(P1, tmax, n, L, [1 2 3]);

% P2 = [0.05 0.05 0.07];
% GM_solve(P2, tmax, n, L, [1 2 3]);

% P3 = [0.4 0.5 0.07];
% GM_solve(P3, tmax, n, L, [1 2 3]);

P4 = [0.05 1.4 0.07];
GM_solve(P4, tmax, n, L, [1 2 3]);

%% Problem 1, Part B
close all
clc
clear

tmax = 100;
n = 100;
L = 10;

P = [0.09 1.25 0.07];
GM_solve(P, tmax, n, L, [1 2 3]);

%% Problem 1, Part C
close all
clc
clear

tmax = 100;
n = 4;
L = 10;

P = [0.05, 1.4, 0.07];

IC_func_sym = @(x, P) IC_func(x, P, n, L);
GM_solve_IC(P, tmax, n, L, [1 2 3], IC_func_sym);


%% Problem 2: 2D pattern formation in the Gray-Scott model
disp('PROBLEM 2')
clear

%% Problem 2, Part A
close all

% rxn_dfsn_gs(0.022, 0.051, 1, 0.5, InitState='square');

% rxn_dfsn_gs(0.01, 0.041, 1, 0.5, InitState='wavefront');

% rxn_dfsn_gs(0.09, 0.059, 1, 0.5, InitState='uSpots');

rxn_dfsn_gs(0.09, 0.059, 1, 0.5, InitState='vSpots');

%% Problem 2, Part B
close all

% rxn_dfsn_gs(0.022, 0.051, 1, 0.7, InitState='square');

rxn_dfsn_gs(0.022, 0.051, 1, 0.8, InitState='square');

%% Problem 2, Part C
close all

% rxn_dfsn_gs(0.01, 0.041, 1, 0.5, InitState='wavefront');

rxn_dfsn_gs(0.01, 0.041, 1, 0.6, InitState='wavefront');


%% Functions

% Problem 1
function y0 = IC_func(x, P, n, L)
    h = ones(size(x));
    a = 1+cos(n*pi*x/L);

    y0 = [a; h];
end
