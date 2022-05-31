close all
clear
clc

addpath('vfit3')
addpath('functions')

% Load dataset
fname='carson';
load([fname '.mat'])

% General parameters
w = 2.*pi.*f; % Frequency (rad/seg)
s = 1j*w; % Complex Frequency
line_length = 1000;
f_Ti=200e3; % Frequency chosen for JMarti model

% Raw impedance and admittance matrices
Z=permute(Ztot_Carson,[3 1 2]);
Y=permute(Ytot_Imag,[3 1 2]);

% Frequency-variant transform matrix from OHLT
    % [Ti,g_dis]=LM_calc_norm_str(ord,freq_siz,Z,Y,f);
    % [Ti,g_dis]=intercheig_QR_decomp(ord,freq_siz,Z,Y);
[Ti,g_dis]=simple_QR_decomp(ord,freq_siz,Z,Y);
Ti=rot_min_imag(Ti,ord,freq_siz);
Ti=list2sqmat(Ti,ord,freq_siz);

% Now let's pick the corresponding matrix for a specific frequency @ f_Ti
interplM2freq = @(M) squeeze(interp1(f,M,f_Ti)); %use interpolation to handle the case where the frequency sample is missing
T=interplM2freq(Ti);
T=real(T); %is this really necessary?
g=interp1(f,g_dis,f_Ti); %grab propagation constant at chosen frequency


