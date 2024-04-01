close all
clc

%% Problem 1: Moran model with no selection (DEMO)
disp('PROBLEM 1')
clear

%% Problem 1, Part A
% initializing parameters
N = 10;
t0 = 1; % day
dt = t0/N;
tf = 50;
tspan = 1:dt:tf;
Nsim = 1000;
n_mutants = zeros(length(tspan), Nsim);

% simulate model
for i = 1:Nsim
    cells = zeros(N,1);
    cells(randi(N)) = 1;
    n_mutants(1, i) = 1;
   
    for j = 2:length(tspan)
        r = randi(N, 2, 1);
        cells(r(1)) = cells(r(2));
        n_mutants(j, i) = sum(cells);
       
        if n_mutants(j, i) == 0
            break
        elseif n_mutants(j, i) == N
            n_mutants(j+1:length(tspan), i) = N;
        end
    end
end

m_mutants = mean(n_mutants, 2);
prob_survival = sum(n_mutants(:,:)' > 0)/Nsim;
prob_survival_tf = prob_survival(end);

figure(1)
hold on
plot(tspan, n_mutants(:, 1:3))
plot(tspan, m_mutants)
xlabel('Time')
ylabel('M')
legend('Sim 1','Sim 2', 'Sim 3', 'Mean of 1000 Sims')
title('PROBLEM 1A: Number of mutant cells M v. time for sims')

figure(2)
hold on
plot(tspan, prob_survival)
xlabel('Time')
ylabel('Probability of survival')
title('PROBLEM 1A: Probability of survival for mutant cancer cells v. time')

fprintf("Probability of mutant survival at final time: %.3f\n", prob_survival_tf)

%% Problem 1, Part B
N = 10;
t0 = 1;
dt = t0/N;
tf = 50;
tspan = 1:dt:tf;
Nsim = 1000;
n_mutants = zeros(length(tspan), Nsim);

for i = 1:Nsim
    n_mutants(1, i) = 1;
    for j = 2:length(tspan)
        prob_increase = n_mutants(j-1, i)/N * (1-n_mutants(j-1, i)/N);
        prob_decrease = n_mutants(j-1, i)/N * (1-n_mutants(j-1, i)/N);
       
        r = rand;
        if r <= prob_increase
            n_mutants(j, i) = n_mutants(j-1, i)+1;
        elseif r<= prob_increase + prob_decrease
            n_mutants(j, i) = n_mutants(j-1, i)-1;
        else
            n_mutants(j, i) = n_mutants(j-1, i);
        end
       
        if n_mutants(j, i) ==0
            break
        elseif n_mutants(j, i) == N
            n_mutants(j+1:length(tspan), i) = N;
            break
        end
    end
end

m_mutants = mean(n_mutants, 2);
prob_survival = sum(n_mutants(:,:)' > 0)/Nsim;
prob_survival_tf = prob_survival(end);

figure(3)
hold on
plot(tspan, n_mutants(:, 1:3))
plot(tspan, m_mutants)
xlabel('Time')
ylabel('M')
legend('Sim 1','Sim 2', 'Sim 3', 'Mean of 1000 Sims')
title('PROBLEM 1B: Number of mutant cells M v. time for sims')

figure(4)
hold on
xlabel('Time')
ylabel('Probability of survival')
plot(tspan, prob_survival)
title('PROBLEM 1B: Probability of survival for mutant cancer cells v. time')

fprintf("Probability of mutant survival at final time: %.3f\n", prob_survival_tf)


%% Problem 2: Moran model with selection
disp('PROBLEM 2')
clear

%% Problem 2, Part A
N = 10;
t0 = 1;
dt = t0/N;
tf = 150;
tspan = 1:dt:tf;
Nsim = 1000;
n_mutants = zeros(length(tspan), Nsim);

s = 0.05;

for i = 1:Nsim
    n_mutants(1, i) = 1;
    for j = 2:length(tspan)
        prob_increase = (1+s) * n_mutants(j-1, i)/N * (1-n_mutants(j-1, i)/N);
        prob_decrease = n_mutants(j-1, i)/N * (1-n_mutants(j-1, i)/N);
       
        r = rand;
        if r <= prob_increase
            n_mutants(j, i) = n_mutants(j-1, i)+1;
        elseif r<= prob_increase + prob_decrease
            n_mutants(j, i) = n_mutants(j-1, i)-1;
        else
            n_mutants(j, i) = n_mutants(j-1, i);
        end
       
        if n_mutants(j, i) ==0
            break
        elseif n_mutants(j, i) == N
            n_mutants(j+1:length(tspan), i) = N;
            break
        end
    end
end

m_mutants = mean(n_mutants, 2);
prob_survival = sum(n_mutants(:,:)' > 0)/Nsim;
prob_survival_tf = prob_survival(end);

figure(5)
hold on
plot(tspan, n_mutants(:, 1:3))
plot(tspan, m_mutants)
xlabel('Time')
ylabel('M')
legend('Sim 1', 'Sim 2', 'Sim 3', 'Mean of 1000 Sims')
title('PROBLEM 2A: Number of mutant cells M v. time for sims')

figure(6)
hold on
xlabel('Time')
ylabel('Probability of survival')
plot(tspan, prob_survival)
title('PROBLEM 2A: Probability of survival for mutant cancer cells v. time')

fprintf("Probability of mutant survival at final time: %.3f\n", prob_survival_tf)

%% Problem 2, Part B
N = 10;
t0 = 1;
g = log10(2)/t0;
s = 0.05;
M0 = 1;
tf = 150;

p = [g, s, N];
[t, M] = ode45(@(t, M) ode2B(t, M, p), [t0 tf], M0);

M_mean = mean(M)
M_150 = M(end)

figure(7)
hold on
plot(t, M, '-b', LineWidth=1.5)
xlabel('Time')
ylabel('M')
title('PROBLEM 2B: Number of mutant cells M v. time using ode45')


%% Problem 2, Part C
N = 100;
t0 = 1;
dt = t0/N;
tf = 150;
tspan = 1:dt:tf;
Nsim = 1000;
n_mutants = zeros(length(tspan), Nsim);

s = 0.05;

for i = 1:Nsim
    n_mutants(1, i) = 1;
    for j = 2:length(tspan)
        prob_increase = (1+s) * n_mutants(j-1, i)/N * (1-n_mutants(j-1, i)/N);
        prob_decrease = n_mutants(j-1, i)/N * (1-n_mutants(j-1, i)/N);
       
        r = rand;
        if r <= prob_increase
            n_mutants(j, i) = n_mutants(j-1, i)+1;
        elseif r<= prob_increase + prob_decrease
            n_mutants(j, i) = n_mutants(j-1, i)-1;
        else
            n_mutants(j, i) = n_mutants(j-1, i);
        end
       
        if n_mutants(j, i) ==0
            break
        elseif n_mutants(j, i) == N
            n_mutants(j+1:length(tspan), i) = N;
            break
        end
    end
end

m_mutants = mean(n_mutants, 2);
prob_survival = sum(n_mutants(:,:)' > 0)/Nsim;
prob_survival_tf = prob_survival(end);

figure(8)
hold on
plot(tspan, n_mutants(:, 1:3))
plot(tspan, m_mutants)
xlabel('Time')
ylabel('M')
legend('Sim 1', 'Sim 2', 'Sim 3', 'Mean of 1000 Sims')
title('PROBLEM 2C: Number of mutant cells M v. time for sims')

figure(9)
hold on
xlabel('Time')
ylabel('Probability of survival')
plot(tspan, prob_survival)
title('PROBLEM 2C: Probability of survival for mutant cancer cells v. time')

fprintf("Probability of mutant survival at final time: %.3f\n", prob_survival_tf)

N = 100;
t0 = 1;
g = log10(2)/t0;
s = 0.05;
M0 = 1;
tf = 150;

p = [g, s, N];
[t, M] = ode45(@(t, M) ode2B(t, M, p), [t0 tf], M0);

M_mean = mean(M)
M_150 = M(end)

figure(10)
hold on
plot(t, M, '-b', LineWidth=1.5)
xlabel('Time')
ylabel('M')
title('PROBLEM 2C: Number of mutant cells M v. time using ode45')


%% Problem 3: Spatial Moran model
disp('PROBLEM 3')
clear

%% Problem 3, Part A
N = 10;
t0 = 1;
dt = t0/N;
tf = 50;
tspan = 1:dt:tf;
Nsim = 1000;

s = 0.05;

n_mutants = zeros(length(tspan), Nsim);

for i = 1:Nsim
    lattice = zeros(1, N);
    lattice(1, randi(N)) = 1;
    n_mutants(1, Nsim) = 1;

    for j = 2:length(tspan)

        idx = randi(N);

        if idx ~= 1 && idx ~= N
            neighbors = [idx-1, idx+1];
            if lattice(neighbors(1)) == 1
                prob = (1+s)/(2+s);
            else
                prob =  1/(2+s);
            end

            if rand < prob
                lattice(idx) = lattice(neighbors(1));
            else
                lattice(idx) = lattice(neighbors(2));
            end

        elseif idx == 1
            lattice(idx) = lattice(idx+1);
        elseif idx == N
            lattice(idx) = lattice(idx-1);
        end
        
        n_mutants(j, i) = sum(lattice);

        if n_mutants(j, i) == 0
            break
        elseif n_mutants(j, i) == N
            n_mutants(j+1:length(tspan), i) = N;
            break
        end
    end
end

m_mutants = mean(n_mutants, 2);
prob_survival = sum(n_mutants(:,:)' > 0)/Nsim;
prob_survival_tf = prob_survival(end);

figure(11)
hold on
plot(tspan, n_mutants(:, 1:3))
plot(tspan, m_mutants)
xlabel('Time')
ylabel('M')
legend('Sim 1', 'Sim 2', 'Sim 3', 'Mean of 1000 Sims')
title('PROBLEM 3A: Number of mutant cells M v. time for sims')

figure(12)
hold on
xlabel('Time')
ylabel('Probability of survival')
plot(tspan, prob_survival)
title('PROBLEM 3A: Probability of survival for mutant cancer cells v. time')

fprintf("Probability of mutant survival at final time: %.3f\n", prob_survival_tf)

%% Problem 3, Part B
N = 100;
t0 = 1;
dt = t0/N;
tf = 50;
tspan = 1:dt:tf;
Nsim = 1000;

s = 0.05;

n_mutants = zeros(length(tspan), Nsim);

for i = 1:Nsim
    lattice = zeros(1, N);
    lattice(1, randi(N)) = 1;
    n_mutants(1, Nsim) = 1;

    for j = 2:length(tspan)

        idx = randi(N);

        if idx ~= 1 && idx ~= N
            neighbors = [idx-1, idx+1];
            if lattice(neighbors(1)) == 1
                prob = (1+s)/(2+s);
            else
                prob =  1/(2+s);
            end

            if rand < prob
                lattice(idx) = lattice(neighbors(1));
            else
                lattice(idx) = lattice(neighbors(2));
            end

        elseif idx == 1
            lattice(idx) = lattice(idx+1);
        elseif idx == N
            lattice(idx) = lattice(idx-1);
        end
        
        n_mutants(j, i) = sum(lattice);

        if n_mutants(j, i) == 0
            break
        elseif n_mutants(j, i) == N
            n_mutants(j+1:length(tspan), i) = N;
            break
        end
    end
end

m_mutants = mean(n_mutants, 2);
prob_survival = sum(n_mutants(:,:)' > 0)/Nsim;
prob_survival_tf = prob_survival(end);

figure(13)
hold on
plot(tspan, n_mutants(:, 1:3))
plot(tspan, m_mutants)
xlabel('Time')
ylabel('M')
legend('Sim 1', 'Sim 2', 'Sim 3', 'Mean of 1000 Sims')
title('PROBLEM 3B: Number of mutant cells M v. time for sims')

figure(14)
hold on
xlabel('Time')
ylabel('Probability of survival')
plot(tspan, prob_survival)
title('PROBLEM 3B: Probability of survival for mutant cancer cells v. time')

fprintf("Probability of mutant survival at final time: %.3f\n", prob_survival_tf)

%% Problem 3, Part C
% No MATLAB code required.

%% Functions

% Problem 2
function dMdt = ode2B(t, M, p)
    g = p(1);
    s = p(2);
    N = p(3);

    g_av = ((M/N)*(1+s)*g) + (((N-M)/N)*g);
    dMdt = (g*M*(1+s)) - (g_av*M);
end
