function [H_mod,F_mod,pol_co] = HF_VF_fun(Ti,g,length,freq,ord)

freq_siz=max(size(freq));

%% Calculation of propagation function
[H_mod,P_mod,H,P,delay]=calc_prop_function(Ti,g,length,ord,freq_siz,freq);

%% Calculation of transfer function
[F_mod,F,f]=calc_transf_function(Ti,g,length,ord,freq);

%% Vector Fitting
[Ti_vf,fit,check]=vector_fit(Ti,ord,freq_siz,freq);
%check

%% Calculation of Poles and Zeros
[pol_co]=calc_poles_zeros(ord,fit);
