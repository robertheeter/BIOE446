close all
clc

%% Problem 1: Simulation of Lotka-Volterra model (DEMO)
disp('PROBLEM 1')
clear

%% Problem 1, Part A
% setting time, species values, and parameters
t = 0;
y = [0.4 0.4];
a = 0.5; b = 3; c = 3; m = 0.1;
p = [a b c m];

dydt = LV_model(t, y, p);
disp(dydt);

%% Problem 1, Part B
% setting timepoints, initial conditions, and parameter values
tspan = linspace(0, 500, 5000);
y0 = [0.4 0.4; 0.3 0.3; 0.2 0.2];

% setting up function for ode simulations
dydt = @(t, y) LV_model(t, y, p);

% initalizing storage for H and P
P = zeros(length(tspan), size(y0, 1));
H = P;

% setting tolerances for ode solver
options = odeset('AbsTol', 1e-8, 'RelTol', 1e-5);

% iterate through the initial conditions and simulate
for ii = 1:size(y0, 1)
    [~, y] = ode15s(dydt, tspan, y0(ii, :), options);
    P(:, ii) = y(:, 2);
    H(:, ii) = y(:, 1);
end

figure(1)
hold on
plot(tspan, H(:,1), LineWidth=1.5)
plot(tspan, P(:,2), LineWidth=1.5)
xlabel('Time')
ylabel('Species abundances')
title('PROBLEM 1B: H, P v. t')
legend('H', 'P', Location='best')
set(gca, 'Yscale', 'log')

figure(2)
hold on
plot(H, P, LineWidth=1.5)
dHdt = @(H, P) a.*H - b.*H.*P;
dPdt = @(H, P) c.*H.*P - m.*P;
fimplicit(dHdt, '--k', LineWidth=1.5);
fimplicit(dPdt, '--k', LineWidth=1.5);
title('PROBLEM 1B: P v. H')
xlabel('H')
ylabel('P')
l = legend('(0.4, 0.4)', '(0.3, 0.3)', '(0.2, 0.2)', location='best');
title(l, '(H_0, P_0)')
set(gca, 'Yscale', 'log')
set(gca, 'Xscale', 'log')

%% Problem 1, Part C
% No MATLAB code required.


%% Problem 2: Introducing a carrying capacity on prey growth
disp('PROBLEM 2')
clear

%% Problem 2, Part A
t = 0;
y = [0.4 0.4];
a = 0.5; b = 3; c = 3; m = 0.1; K = 0.4;
p = [a b c m K];

dydt = carrying_capacity(t, y, p);
disp(dydt);

%% Problem 2, Part B
tspan = linspace(0, 500, 5000);
y0 = [0.4 0.4; 0.3 0.3; 0.2 0.2];

dydt = @(t, y) carrying_capacity(t, y, p);

P = zeros(length(tspan), size(y0, 1));
H = P;

options = odeset('AbsTol', 1e-8, 'RelTol', 1e-5);

for ii = 1:size(y0, 1)
    [~, y] = ode15s(dydt, tspan, y0(ii, :), options);
    P(:, ii) = y(:, 2);
    H(:, ii) = y(:, 1);
end

figure(3)
hold on
plot(tspan, H(:,1), LineWidth=1.5)
plot(tspan, P(:,2), LineWidth=1.5)
xlabel('Time')
ylabel('Species abundances')
title('PROBLEM 2B: H, P v. t')
legend('H', 'P', Location='best')
set(gca, 'Yscale', 'log')

figure(4)
hold on
plot(H, P, LineWidth=1.5)
dHdt = @(H, P) a.*H.*(1-(H./K)) - b.*H.*P;
dPdt = @(H, P) c.*H.*P - m.*P;
fimplicit(dHdt, '--k', LineWidth=1.5);
fimplicit(dPdt, '--k', LineWidth=1.5);
title('PROBLEM 2B: P v. H')
xlabel('H')
ylabel('P')
l = legend('(0.4, 0.4)', '(0.3, 0.3)', '(0.2, 0.2)', location='best');
title(l, '(H_0, P_0)')
set(gca, 'Yscale', 'log')
set(gca, 'Xscale', 'log')


%% Problem 3: Simulating a Macarthur-Rozenzweig model
disp('PROBLEM 2')
clear

%% Problem 3, Part A
t = 0;
y = [0.4 0.4];
a = 0.5; b = 3; eps = 1; tau = 3; m = 0.1; K = 0.4;
p = [a b eps tau m K];

dydt = MR_model(t, y, p);
disp(dydt);

%% Problem 3, Part B
tspan = linspace(0, 500, 5000);
y0 = [0.4 0.4; 0.3 0.3; 0.2 0.2];

dydt = @(t, y) MR_model(t, y, p);

P = zeros(length(tspan), size(y0, 1));
H = P;

options = odeset('AbsTol', 1e-8, 'RelTol', 1e-5);

for ii = 1:size(y0, 1)
    [~, y] = ode15s(dydt, tspan, y0(ii, :), options);
    P(:, ii) = y(:, 2);
    H(:, ii) = y(:, 1);
end

figure(5)
hold on
plot(tspan, H(:,1), LineWidth=1.5)
plot(tspan, P(:,2), LineWidth=1.5)
xlabel('Time')
ylabel('Species abundances')
title('PROBLEM 3B: H, P v. t')
legend('H', 'P', Location='best')
set(gca, 'Yscale', 'log')

figure(6)
hold on
plot(H, P, LineWidth=1.5)
dHdt = @(H, P) a.*H.*(1-(H./K)) - (b.*H.*P)./(1+(b.*tau.*H));
dPdt = @(H, P) (eps.*b.*H.*P)./(1+(b.*tau.*H)) - m.*P;
fimplicit(dHdt, '--k', LineWidth=1.5);
fimplicit(dPdt, '--k', LineWidth=1.5);
title('PROBLEM 3B: P v. H')
xlabel('H')
ylabel('P')
l = legend('(0.4, 0.4)', '(0.3, 0.3)', '(0.2, 0.2)', location='best');
title(l, '(H_0, P_0)')
set(gca, 'Yscale', 'log')
set(gca, 'Xscale', 'log')


%% Functions

% Problem 1
function dydt = LV_model(t, y, p)
% Rate laws for standard LV model
    % setting species
    H = y(1);
    P = y(2);
    
    % setting parameters
    a = p(1); b = p(2); c = p(3); m = p(4);
    
    % calculating rates of change
    dydt = zeros(size(y));
    
    % dHdt
    dydt(1) = a*H - b*H*P;
    
    % dPdt
    dydt(2) = c*H*P - m*P;
end

% Problem 2
function dydt = carrying_capacity(t, y, p)
    H = y(1);
    P = y(2);

    a = p(1); b = p(2); c = p(3); m = p(4); K = p(5);

    dydt = zeros(size(y));
    dydt(1) = a*H*(1-(H/K)) - b*H*P;
    dydt(2) = c*H*P - m*P;
end

% Problem 3
function dydt = MR_model(t, y, p)
    H = y(1);
    P = y(2);

    a = p(1); b = p(2); eps = p(3); tau = p(4); m = p(5); K = p(6);
    
    dydt = zeros(size(y));
    dydt(1) = a*H*(1-(H/K)) - (b*H*P)/(1+(b*tau*H));
    dydt(2) = (eps*b*H*P)/(1+(b*tau*H)) - m*P;
end
