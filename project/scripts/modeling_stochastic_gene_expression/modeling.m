%% MODELING STOCHASTIC GENE EXPRESSION

%% EXPERIMENTAL DATA
total_cells = 196; % total number of cells, from image analysis
total_mRNA_f = 93885.9915; % total free mRNA
total_mRNA_n = 1405083.355; % total nascent mRNA

%% CALCULATING VARIABLES
mRNA_f = total_mRNA_f/total_cells % mean free mRNA
mRNA_n = total_mRNA_n/total_cells % mean nascent mRNA

gamma_1 = 0.05 % degradation rate, min^-1, given
k_tsc = 21 % transcription rate, min^-1, given

k_release = (mRNA_f/mRNA_n)*gamma_1 % release rate, min^-1

k_prod = mRNA_n*k_release % production rate, min^-1

ratio = (k_prod/k_tsc)/(1-(k_prod/k_tsc)); % ratio k_on/k_off

%% FINDING K_ON AND K_OFF
k_on = 0.05 % arbitrary, min^-1
k_off = k_on/ratio % min^-1

project = sbioloadproject('GeneExpression.sbproj');
model = project.m1;
model.Parameters(1).Value = k_on;
model.Parameters(2).Value = k_off;
model.Parameters(4).Value = k_release;

[success, out, modek] = sbiosteadystate(model, 'method', 'simulation');
model.Species(1).Value = 1;
model.Species(2).Value = 0;

config = getconfigset(model);
config.StopTime = 100;
config.SolverType = 'ssa';

ens_data = sbioensemblerun(model, 1000);

mRNA_n_sim = zeros(1, length(ens_data));
mRNA_f_sim = zeros(1, length(ens_data));

for i = 1:length(ens_data)
    mRNA_n_sim(i) = ens_data(i).Data(end, 3);
    mRNA_f_sim(i) = ens_data(i).Data(end, 4);
end

figure(1)
hold on
histogram(mRNA_n_sim)
xlabel('Nascent mRNA [#] at t=100')
ylabel('Frequency')
title('Histogram of nascent mRNA at t_{max} = 100 min')

figure(2)
hold on
histogram(mRNA_f_sim)
xlabel('Free mRNA [#] at t=100')
ylabel('Frequency')
title('Histogram of free mRNA at t_{max} = 100 min')

%% DETERMINE STATISTICS
mRNA_n_sim_var = var(mRNA_n_sim)

total_mRNA_sim = mRNA_n_sim + mRNA_f_sim;

total_mRNA_sim_mean = mean(total_mRNA_sim)
total_mRNA_sim_var = var(total_mRNA_sim)
total_mRNA_sim_cv = sqrt(total_mRNA_sim_var)/total_mRNA_sim_mean
