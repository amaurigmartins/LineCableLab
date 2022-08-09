% Time Domain Calculation
function v=tm_dmn_clc(V,f,e,data_t_sim,c,ord)

Pn=exp(c*data_t_sim)*e;

%sigma=1; % No damping window, FFT
%sigma=0.5*(1+cos(pi*f/(e/2))); %Hanning
%sigma=sinc(f/(e/2)); %Lanczos
sigma=0.42+0.5*cos(pi*f/(e/2))+0.08*cos(2*pi*f/(e/2)); %Blackman
%sigma=1-(f/(e/2)).^2; %Riesz

for o=1:1:2*ord
    V(o,:)=V(o,:).*sigma;
end

hsmp=max(size(V));

for p=1:1:hsmp
    V(:,hsmp*2-p)=V(:,p); % Δημιουργία του ειδώλου στη τάση V
end

v=ifft(transpose(V),'symmetric');% IFFT θεωρώντας τη τάση συζυγή συμμετρική - (2*max(size(f)) x 2*ord)

for o=1:1:2*ord
    v(:,o)=real(transpose(Pn).*v(:,o));
end