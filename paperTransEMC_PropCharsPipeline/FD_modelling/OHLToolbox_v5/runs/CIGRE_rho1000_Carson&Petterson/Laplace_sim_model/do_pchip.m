function [gamma,Z,Y,Ti,Ys,Yr]=do_pchip(ord,freq,f,gamma_dis,Z_dis,Y_dis,Ti_dis,Ys_dis,Yr_dis)

% gamma_real=zeros(max(size(f)),ord);
% gamma_imag=zeros(max(size(f)),ord);
gamma=zeros(max(size(f)),ord);

% Ys_real=zeros(max(size(f)),ord);
% Ys_imag=zeros(max(size(f)),ord);
Ys=zeros(max(size(f)),ord);

% Yr_real=zeros(max(size(f)),ord);
% Yr_imag=zeros(max(size(f)),ord);
Yr=zeros(max(size(f)),ord);

% Z_real=zeros(max(size(f)),ord^2);
% Z_imag=zeros(max(size(f)),ord^2);
Z=zeros(max(size(f)),ord^2);

% Y_real=zeros(max(size(f)),ord^2);
% Y_imag=zeros(max(size(f)),ord^2);
Y=zeros(max(size(f)),ord^2);

% Ti_real=zeros(max(size(f)),ord^2);
% Ti_imag=zeros(max(size(f)),ord^2);
Ti=zeros(max(size(f)),ord^2);

for o=1:1:ord
    gamma_real(:,o)=pchip(freq,real(gamma_dis(:,o)));
    gamma_imag(:,o)=pchip(freq,imag(gamma_dis(:,o)));
    
    Ys_real(:,o)=pchip(freq,real(Ys_dis(:,o)));
    Ys_imag(:,o)=pchip(freq,imag(Ys_dis(:,o)));
    
    Yr_real(:,o)=pchip(freq,real(Yr_dis(:,o)));
    Yr_imag(:,o)=pchip(freq,imag(Yr_dis(:,o)));
end

for o=1:1:ord^2
    Z_real(:,o)=pchip(freq,real(Z_dis(:,o)));
    Z_imag(:,o)=pchip(freq,imag(Z_dis(:,o)));
    
    Y_real(:,o)=pchip(freq,real(Y_dis(:,o)));
    Y_imag(:,o)=pchip(freq,imag(Y_dis(:,o)));
    
    Ti_real(:,o)=pchip(freq,real(Ti_dis(:,o)));
    Ti_imag(:,o)=pchip(freq,imag(Ti_dis(:,o)));
end

for o=1:1:ord
    gamma(:,o)=ppval(gamma_real(:,o),f)+1j*ppval(gamma_imag(:,o),f);
    Ys(:,o)=ppval(Ys_real(:,o),f)+1j*ppval(Ys_imag(:,o),f);
    Yr(:,o)=ppval(Yr_real(:,o),f)+1j*ppval(Yr_imag(:,o),f);
end

for o=1:1:ord^2
    Z(:,o)=ppval(Z_real(:,o),f)+1j*ppval(Z_imag(:,o),f);
    Y(:,o)=ppval(Y_real(:,o),f)+1j*ppval(Y_imag(:,o),f);
    Ti(:,o)=ppval(Ti_real(:,o),f)+1j*ppval(Ti_imag(:,o),f);
end
