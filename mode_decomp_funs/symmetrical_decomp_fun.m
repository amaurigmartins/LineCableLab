function [outputArg1,outputArg2] = symmetrical_decomp_fun(matriz z cheia,matriz y cheia, vetor de frequencias)
%SYMMETRICAL_DECOMP_FUN Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;

% atenção para redução de kron - nao sabemos a priori quais cabos são CPR
cpr_idx=[7 8];

%fortescue
A=exp(1i*2*pi/3);
F=[1 1 1; 1 A^2 A; 1 A A^2];
mZ012=(F\Z)*F;

% clarke 
% F=(1/3)*[2 -1 -1; 0 sqrt(3) -sqrt(3); 1 1 1];
% T=inv(F);
% mZ012=inv(T)*Z*T;

%retornar matrizes Z e Y modais e parametros de propagação

for k=1:n_freq_samples

end

end

