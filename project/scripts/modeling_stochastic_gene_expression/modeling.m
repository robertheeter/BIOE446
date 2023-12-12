% MODELING STOCHASTIC GENE EXPRESSION

%% PARAMETERS
y_1 = 0.05; % degradation rate, min^-1, given
k_tsc = 21; % transcription rate, min^-1, given

mRNA_f = 0; % mean free mRNA, from image analysis
mRNA_n = 0; % mean nascent mRNA, from image analysis

k_release = (mRNA_f/mRNA_n)*y_1; % release rate, min^-1

k_prod = mRNA_n*k_release; % production rate, min^-1

ratio = 1/((k_tsc/k_prod)-1); % ratio k_on/k_off

k_on = 1; % arbitrary, min^-1
k_off = k_on/ratio; % min^-1
