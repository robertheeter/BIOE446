close all
clc

%% Challenge Problem 1
disp('CHALLENGE PROBLEM 1')
clear

x = rand;
Tnext = (1/5)*log(1/x);
Rind = rand;

if Rind < 2/5
    n = 4;
else
    n = 2;
end
disp(n)


%% Problem 1: Simulation of stochastic gene expression dynamics (DEMO)
disp('PROBLEM 1')
clear

%% Problem 1, Part A
% initialize
M0 = 0;
alpha = 30;
gamma = 1;
tmax = 720;

i = 1;
M(1) = 0;
t(1) = 0;

t_current = 0;
m_current = 0;

% SSA
while t_current < tmax
    % calculate production and degradation
    production = alpha;
    degradation = gamma*m_current;
    total_rate = production + degradation;

    % determine wait time until the next event
    tau = exprnd(1/total_rate);
    t_current = t_current + tau;

    % determine which event occurred (i.e. birth or death)
    if (rand < production/total_rate)
        m_current = m_current + 1;
    else
        m_current = m_current - 1;
    end

    i = i+1;
    M(i) = m_current;
    t(i) = t_current;
end

% ode solution
[t_ode, M_ode] = ode45(@(t,m) alpha-gamma*m, [0 10], 0);

% plot results
figure(1)
hold on
plot(t(t<=10), M(t<=10), LineWidth=1.5, DisplayName='SSA')
plot(t_ode, M_ode, LineWidth=1.5, DisplayName='ODE45')
xlabel('Time [min]')
ylabel('mRNA [#]')
title('Problem 1A')
legend(location='best')
legend box off

%% Problem 1, Part B
M_reg = interp1(t, M, 0:0.1:720, 'previous');
mean_M_reg = mean(M_reg)
var_M_reg = var(M_reg)


%% Problem 2: Modeling stochastic gene expression with SimBiology
disp('PROBLEM 2')
clear

%% Problem 2, Part A
simbio = load('problem2a.mat');

figure(2)
hold on
plot(simbio.problem2a.Time, simbio.problem2a.Data, LineWidth=1.5, DisplayName='Stochastic')
xlabel('Time [min]')
ylabel('mRNA [#]')
title('Problem 2A')
legend(location='best')
legend box off

%% Problem 2, Part B
t = simbio.problem2a.Time;
M = simbio.problem2a.Data;
M_reg = interp1(t, M, 0:0.1:720, 'previous');
mean_M_reg = mean(M_reg)
var_M_reg = var(M_reg)

%% Problem 2, Part C
poisson_dist = poisspdf(0:720, mean(M_reg));

% plot histogram for M
figure(3)
hold on
histogram(M_reg, 'normalization', 'pdf', DisplayName='Simulation')
plot(poisson_dist, LineWidth=1.5, DisplayName='Poisson')
xlim([0 50])
xlabel('mRNA [#]')
ylabel('Frequency')
title('Problem 2C')
legend(location='best')

%% Problem 2, Part D
project = sbioloadproject('Birth-Death.sbproj');
model = project.m1;
config = getconfigset(model);
set(config, 'SolverType', 'ssa');

M_ens = sbioensemblerun(model, 1000);

y = zeros(1000, 1);
for i = 1:1000
    y(i) = M_ens(i).Data(end);
end

% plot results
figure(4)
hold on
histogram(y)
xlabel('mRNA [#] at t=720')
ylabel('Frequency')
title('Problem 2D: Histogram of M at t_{max} = 720 min')


%% Problem 3: Expanding the transcriptional bursting model
disp('PROBLEM 3')
clear

%% Problem 3, Part A
simbio = load('problem3a.mat');
X = simbio.problem3a.Data(:,4);
mean_X = mean(X)

%% Problem 3, Part B
project = sbioloadproject('Multiple-Promoter.sbproj');
model = project.m1;
config = getconfigset(model);
set(config, 'SolverType', 'ssa');

X_ens = sbioensemblerun(model, 1000);

%% Problem 3, Part C
X = zeros(1000, 1);
for i = 1:1000
    data = X_ens(i).Data(:,4);
    X(i) = data(end);
end

% plot results
figure(5)
hold on
histogram(X)
xlabel('Protein X [#] at t=240')
ylabel('Frequency')
title('Problem 3C: Histogram of X at t_{max} = 240 min')

%% Problem 3, Part D
mean_X = mean(X)
stdev_X = std(X)
co_var_X = stdev_X/mean_X

%% Problem 3, Part E
simbio = load('problem3e.mat');
X = simbio.problem3e.Data(:,4);
mean_X = mean(X)

project = sbioloadproject('Multiple-Promoter.sbproj');
model = project.m1;
config = getconfigset(model);
set(config, 'SolverType', 'ssa');

X_ens = sbioensemblerun(model, 1000);

X = zeros(1000, 1);
for i = 1:1000
    data = X_ens(i).Data(:,4);
    X(i) = data(end);
end

% plot results
figure(6)
hold on
histogram(X)
xlabel('Protein X [#] at t=240')
ylabel('Frequency')
title('Problem 3E: Histogram of X at t_{max} = 240 min')


%% Problem 4: Feedback and regulatory noise
disp('PROBLEM 4')
clear

%% Problem 4, Part A
% No MATLAB code required.

%% Problem 4, Part B
project = sbioloadproject('Feedback.sbproj');
model = project.m1;
config = getconfigset(model);
set(config, 'SolverType', 'ssa');

X_ens = sbioensemblerun(model, 1000);

%% Problem 4, Part C
X = zeros(1000, 1);
for i = 1:1000
    data = X_ens(i).Data(:,4);
    X(i) = data(end);
end

% plot results
figure(7)
hold on
histogram(X)
xlabel('Protein X [#] at t=240')
ylabel('Frequency')
title('Problem 4C: Histogram of X at t_{max} = 240 min (positive feedback)')

mean_X = mean(X)
stdev_X = std(X)
co_var_X = stdev_X/mean_X

%% Problem 4, Part D
project = sbioloadproject('Feedback.sbproj');
model = project.m1;
config = getconfigset(model);
set(config, 'SolverType', 'ssa');

X_ens = sbioensemblerun(model, 1000);

%% Problem 4, Part E
X = zeros(1000, 1);
for i = 1:1000
    data = X_ens(i).Data(:,4);
    X(i) = data(end);
end

% plot results
figure(8)
hold on
histogram(X)
xlabel('Protein X [#] at t=240')
ylabel('Frequency')
title('Problem 4E: Histogram of X at t_{max} = 240 min (negative feedback)')

mean_X = mean(X)
stdev_X = std(X)
co_var_X = stdev_X/mean_X


%% Functions

% None
