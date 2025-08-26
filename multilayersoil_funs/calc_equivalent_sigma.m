function [ sigma_eq ] = calc_equivalent_sigma( sigma, h, freq )

Nmax=length(sigma);
sigma_eq = 0;
mu0=4*pi*1e-7;

f=@(sigma_top,sigma_bottom,h) sigma_top * (((sqrt(sigma_top)+sqrt(sigma_bottom))-(sqrt(sigma_top)-sqrt(sigma_bottom))*exp(-2*h*sqrt(pi*freq*mu0*sigma_top)))/ ...
    ((sqrt(sigma_top)+sqrt(sigma_bottom))+(sqrt(sigma_top)-sqrt(sigma_bottom))*exp(-2*h*sqrt(pi*freq*mu0*sigma_top))))^2;

for n=Nmax-1:-1:1
    if n==Nmax-1
        sigma_eq = f(sigma(n),sigma(n+1),h(n));
    else
        sigma_eq = f(sigma(n),sigma_eq,h(n));
    end
end

fprintf('Equivalent conductivity (sigma_eq): %g S/m\n', sigma_eq);

end