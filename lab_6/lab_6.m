close all
clc

%% Problem 1: Simulation of a single generation of cell division and mutation (DEMO)
disp('PROBLEM 1')
clear

%% Problem 1, Part A
% defining variables
N = 100;
mu = 1e-3;
N_mutants = 0;

% calculating expected number of mutants
expected_mutants = N*mu;
fprintf('PROBLEM 1A: The expected number of mutants is %i.\n', expected_mutants);

% simulating whether or not mutation occurs during replication
for ii = 1:N
    if rand < mu
        N_mutants = N_mutants + 1;
    end
end

fprintf('PROBLEM 1A: The simulated number of mutants is %i.\n', N_mutants)

%% Problem 1, Part B
% defining variables
C = 100;
N = 2000;
mu = 1e-3;
N_mutants = zeros(1, C);

% iterating through each colony
for ii = 1:C
    % iterating through each cell to see if a mutation occurs during
    % division
    for jj = 1:N
        if rand < mu
            N_mutants(ii) = N_mutants(ii) + 1;
        end
    end
end

figure(1)
histogram(N_mutants);
title('Distribution of number of mutants')
xlabel('Number of mutants')
ylabel('Frequency')

avg = mean(N_mutants)
variance = var(N_mutants)

%% Problem 1, Part C
% defining variables
C = 4000;
N = 2000;
mu = 1e-3;
N_mutants = zeros(1, C);

% random probability from Poisson distribution
for ii = 1:C
    n_mutations = poissrnd(N*mu);
    N_mutants(ii) = n_mutations;
end

figure(2)
histogram(N_mutants);
title('Distribution of number of mutants using Poisson distribution')
xlabel('Number of mutants')
ylabel('Frequency')

avg = mean(N_mutants)
variance = var(N_mutants)


%% Problem 2: Simulation of many generations of cell division and mutation
disp('PROBLEM 2')
clear

%% Problem 2, Part A
g = 15;
N = 500;
mu = 1e-7;
N_mutants = zeros(1, g+1);

for ii = 2:g+1
    n_mutations = poissrnd(N*mu);
    N_mutants(ii) = N_mutants(ii-1)*2 + n_mutations;
    N = N*2;
end

disp(N_mutants)

figure(3)
plot(0:g, N_mutants, LineWidth=1.5)
title('Number of mutants in each generation')
xlabel('Generation number')
ylabel('Number of mutants')

%% Problem 2, Part B
g = 15;
C = 1000;
mu = 1e-7;
N_mutants = zeros(C, g+1);

for ii = 1:C
    N = 500;
    for jj = 2:g+1
        n_mutations = poissrnd(N*mu);
        N_mutants(ii, jj) = N_mutants(ii, jj-1)*2 + n_mutations;
        N = N*2;
    end
end

N_mutants_end = N_mutants(:, end);

figure(4)
histogram(N_mutants_end);
title('Distribution of number of mutants')
xlabel('Number of mutants')
ylabel('Frequency')

avg = mean(N_mutants_end)
variance = var(N_mutants_end)

%% Problem 2, Part C
[B, I] = maxk(N_mutants_end, 5);
C_largest_N_mutants = N_mutants(I, :);

N_mutants_other = N_mutants;
N_mutants_other(I, :) = [];
C_random = N_mutants_other(randi([1 (C-5)], 1, 5), :);

figure(5)
hold on
plot(0:g, C_random, '-b', LineWidth=1.5, DisplayName='C_{random}')
plot(0:g, C_largest_N_mutants, '-m', LineWidth=1.5, DisplayName = 'C_{largest}')
title('Colonies N_{mutants}')
xlabel('Generation number')
ylabel('Number of mutants')
legend(Location="best")
legend box off

%% Problem 3: Inference of mutation rate from experimental data
disp('PROBLEM 3')
% clear

%% Problem 3, Part A
[f,~] = size(find(N_mutants(1:100,15)==0));
f0 = f/100;
u = -log(f0)/N

%% Problem 3, Part B
C = 100;
m = mean(N_mutants(1:100,15));
fsolve(@(mu) mu*N*log(C*mu*N) - m, 1e-7)


%% Functions

% None
