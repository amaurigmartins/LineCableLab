clear all
clc
a=1;

Zsolo = v2_calc_impedancia_solo_v5_gamma_equiv(a);

% Funcao onde entran os dados dos condutores
function [xcondut,ycondut,ncondut,raio,RMG,mizero,sigma,resist,epslonzero,ee,Rdc,freq,Nbund]=v2_dados_condutores(v_fic)
%% Line Geometry
% column 1—conductor number
% column 2-- x position of each cond in m
% column 3-- y position of each cod - Hmax in m
% column 4-- y position of each cod  -Hmin in m
% column 5-- radii of each conductor
% column 6-- number of conductor in bundle
% column 7-- distance between conductors in bundle
% column 8-- conductor resistivity - Resistivity in Ohm-m (no es necesario)
% column 9-- Resistencia DC [ohm/km]

Geom=[1  -2.9 25.32 20.72 0.0160  1 0.4 5.0668e-08 0.063% 1 CONDUCTOR POR BUNDLE: CASO 1
      2   2.9 23.46 19.86 0.0160  1 0.4 5.0668e-08 0.063
      3  -2.9 21.60 16.00 0.0160  1 0.4 5.0668e-08 0.063
      4   0.0 30.00 27.00 0.0079  1 0.4 7.2545e-08 0.500
      ];


Geom(:,3)=(1/3)*Geom(:,3)+(2/3)*Geom(:,4);
ncondut = Geom(max(Geom(:,1)),1);
xcondut=Geom(:,2).';
ycondut=Geom(:,3).';
raio=Geom(:,5);
Nbund=Geom(:,6);
espac=Geom(:,7);
k4=espac./(2.*sin(pi./Nbund));
RMG=(Nbund.*raio.*(k4.^(Nbund-1))).^(1./Nbund);

Rdc=Geom(:,9);
sigma=1./(pi*(raio.^2).*Rdc/1e3);
%%
% permeabilidade  do meio onde a linha está inserida [H/m]
mizero = 4*pi*1e-7; 
% permissividade dielétrica do meio onde a linha está inserida [F/km]
epslonzero= 8.854187817600001e-12; 
% faixa de frequencia
ee=-2:1/30:6;
freq=10.^ee;
% Resistividade do solo [Ohms x m]
resist = 10000;
end



% Funcao onde calculo a Z do solo considerando multicamadas e formulação do
% Nakagawa
function [zsolo] = v2_calc_impedancia_solo_v5_gamma_equiv(vvv)

v_fic=1;
[xcondut,ycondut,ncondut,~,~,mizero,~,resist,epslonzero,~,~,freq]=v2_dados_condutores(v_fic);
%Coordenadas dos condutores imagens
for j = 1: ncondut
   xi(j) =  xcondut(j);
   yi(j) = -ycondut(j);
end

%% Dados das camadas
dd=[2;4;6]; %% altura das camadas

epslonr=[10;10;10;10];

p0=[resist;5000;1000;500];


%%
camadas=size(p0,1);

% Solo com Parametros constantes
sigma_g=1./p0;
epslonr=epslonr.*ones(camadas,length(freq));
sigma_g=sigma_g.*ones(camadas,length(freq));
        
%% Impedancia solo
p1 = 1i*2*freq*mizero;

gamma0_quad = -mizero*epslonzero*(2*pi*freq).^2;
gamma1_quad = 2i*pi*freq*mizero.*(2i*pi*freq.*epslonzero.*epslonr+sigma_g);
          
gamma1_e_matrix=sqrt(gamma1_quad-gamma0_quad);
 

ii=size(gamma1_e_matrix,1);
while 1
    
    if ii==size(gamma1_e_matrix,1)
        
        num=(gamma1_e_matrix(ii-1,:)+gamma1_e_matrix(ii,:)-(gamma1_e_matrix(ii-1,:)-gamma1_e_matrix(ii,:)).*exp(-2*dd(ii-1)*gamma1_e_matrix(ii-1,:)));
        den=(gamma1_e_matrix(ii-1,:)+gamma1_e_matrix(ii,:)+(gamma1_e_matrix(ii-1,:)-gamma1_e_matrix(ii,:)).*exp(-2*dd(ii-1)*gamma1_e_matrix(ii-1,:)));
        gamma1_equiv_matrix=gamma1_e_matrix(ii-1,:).*(num./den);
      
        if ii==2
            break
        end
        
    elseif ii~=[size(gamma1_e_matrix,1),2]
        num=(gamma1_e_matrix(ii-1,:)+gamma1_equiv_matrix-(gamma1_e_matrix(ii-1,:)-gamma1_equiv_matrix).*exp(-2*dd(ii-1)*gamma1_e_matrix(ii-1,:)));
        den=(gamma1_e_matrix(ii-1,:)+gamma1_equiv_matrix+(gamma1_e_matrix(ii-1,:)-gamma1_equiv_matrix).*exp(-2*dd(ii-1)*gamma1_e_matrix(ii-1,:)));
        gamma1_equiv_matrix=gamma1_e_matrix(ii-1,:).*(num./den);
    elseif ii==2
        num=(gamma1_e_matrix(ii-1,:)+gamma1_equiv_matrix-(gamma1_e_matrix(ii-1,:)-gamma1_equiv_matrix).*exp(-2*dd(ii-1)*gamma1_e_matrix(ii-1,:)));
        den=(gamma1_e_matrix(ii-1,:)+gamma1_equiv_matrix+(gamma1_e_matrix(ii-1,:)-gamma1_equiv_matrix).*exp(-2*dd(ii-1)*gamma1_e_matrix(ii-1,:)));
        gamma1_equiv_matrix=gamma1_e_matrix(ii-1,:).*(num./den);
        break
    end
ii=ii-1;    
    
end

gamma_equiv_quadra=gamma1_equiv_matrix.^2;

%% S1
f1 = @(x,j) exp(-2*ycondut(j)*x);
f3 = @(x) (x + sqrt(x.^2 + gamma_equiv_quadra));

%%
ph=1;
for j = 1:ncondut
    for k = 1:ncondut
        if j==k % Zsolo_proprio
            fw = @(x) f1(x,j)./f3(x);
            zsolo_S1 = p1.*integral(fw, 0, Inf, 'arrayvalued', true)*1000; %%S1
 
            S1_int(:,ph)=zsolo_S1.';
        ph=ph+1;
        else % Zsolo_mutuo
        dx = xcondut(j) - xi(k);
        f4 = @(x,j,k,dx) (exp(-(ycondut(j)+ycondut(k))*x)).*cos(dx*x);
        fw = @(x) f4(x,j,k,dx)./f3(x);
        zsolo_S1 = p1.*integral(fw, 0, Inf, 'arrayvalued', true)*1000; %%S1

        S1_int(:,ph)=zsolo_S1.';       
        ph=ph+1;
        end
    end
end
%%
      
zsolo_v2=S1_int;

zsolo = reshape(zsolo_v2.',ncondut,ncondut,[]);
   
end




