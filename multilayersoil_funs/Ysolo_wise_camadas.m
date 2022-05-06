clear all
clc
a=1;

Ysolo = v2_calc_admittance_soil_v5_gamma_equiv(a);

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



% Funcao onde calculo a Y do solo considerando multicamadas e formulação do
% Nakagawa
function [Y_freq]=v2_calc_admittance_soil_v5_gamma_equiv(v_fic)

v_fic=1;
[xcondut,ycondut,ncondut,raio,RMG,mizero,sigma,resist,epslonzero,~,~,freq]=v2_dados_condutores(v_fic);
%Coordenadas dos condutores imagens
for j = 1: ncondut
   xi(j) =  xcondut(j);
   yi(j) = -ycondut(j);
end

MODELO_Y=2;1; %#1->só Externo, #2->Externo+Solo (aproximado_Vi), #3->Externo+Solo (Completo_Vii)

%%
%matriz de coeficientes de potencial
for j = 1:ncondut
   for k = j:ncondut
      
      if j == k
         dezinho = RMG(j);

      else
         dx = xcondut(j) - xcondut(k);
         dy = ycondut(j) - ycondut(k);
         dezinho = sqrt(dx^2 + dy^2);
      end
      
      dx = xcondut(j) - xi(k);
      dy = ycondut(j) - yi(k);
      dezao = sqrt(dx^2 + dy^2);
      
      pot(j,k) = log(dezao/dezinho);
      pot(k,j) = pot(j,k);
   end
end

pot0=pot*1e-3;

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
   
%% Admitancia solo
p1 = 2*ones(1,length(freq));

gamma1_quad = 2i*pi*freq*mizero.*(2i*pi*freq.*epslonzero.*epslonr+sigma_g);  
gamma0_quad = -mizero*epslonzero*(2*pi*freq).^2;
gamma1_e_matrix=sqrt(gamma1_quad-gamma0_quad);
%%
ii=size(gamma1_e_matrix,1);
while 1
    
    if ii==size(gamma1_e_matrix,1)
        
        num=(gamma1_e_matrix(ii-1,:)+gamma1_e_matrix(ii,:)-(gamma1_e_matrix(ii-1,:)-gamma1_e_matrix(ii,:)).*exp(-2*dd(ii-1)*gamma1_e_matrix(ii-1,:)));
        den=(gamma1_e_matrix(ii-1,:)+gamma1_e_matrix(ii,:)+(gamma1_e_matrix(ii-1,:)-gamma1_e_matrix(ii,:)).*exp(-2*dd(ii-1)*gamma1_e_matrix(ii-1,:)));
        gamma1_equiv_matrix=gamma1_e_matrix(ii-1,:).*(num./den);
      
        if ii==2
            break
        end
        
    elseif ii~= [size(gamma1_e_matrix,1),2]
        num=(gamma1_e_matrix(ii-1,:)+gamma1_equiv_matrix-(gamma1_e_matrix(ii-1,:)-gamma1_equiv_matrix).*exp(-2*dd(ii-1)*gamma1_e_matrix(ii-1,:)));
        den=(gamma1_e_matrix(ii-1,:)+gamma1_equiv_matrix+(gamma1_e_matrix(ii-1,:)-gamma1_equiv_matrix).*exp(-2*dd(ii-1)*gamma1_e_matrix(ii-1,:)));
        gamma1_equiv_matrix=gamma1_e_matrix(ii-1,:).*(num./den);
    elseif ii==[2]
        num=(gamma1_e_matrix(ii-1,:)+gamma1_equiv_matrix-(gamma1_e_matrix(ii-1,:)-gamma1_equiv_matrix).*exp(-2*dd(ii-1)*gamma1_e_matrix(ii-1,:)));
        den=(gamma1_e_matrix(ii-1,:)+gamma1_equiv_matrix+(gamma1_e_matrix(ii-1,:)-gamma1_equiv_matrix).*exp(-2*dd(ii-1)*gamma1_e_matrix(ii-1,:)));
        gamma1_equiv_matrix=gamma1_e_matrix(ii-1,:).*(num./den);
        break
    end
ii=ii-1;    
    
end

gamma_equiv_quadra=gamma1_equiv_matrix.^2;

n_gamma_equiv_quadra=gamma_equiv_quadra./gamma0_quad+1;

%% S2
f1 = @(x,j) exp(-2*ycondut(j)*x);
f5 = @(x)  (x.*n_gamma_equiv_quadra + sqrt(x.^2 + gamma_equiv_quadra));


ph=1;
for j = 1:ncondut
    for k = 1:ncondut
        if j==k % Ysolo_proprio
            fw2 = @(x) f1(x,j)./f5(x);
            pgg_S2 = p1.*integral(fw2, 0, Inf, 'arrayvalued', true);
   
            pg_S2(:,ph)=pgg_S2.';
        ph=ph+1;
        else % Ysolo_mutuo
        dx = xcondut(j) - xi(k);
        f4 = @(x,j,k,dx) (exp(-(ycondut(j)+ycondut(k))*x)).*cos(dx*x);
       
        fw2 = @(x) f4(x,j,k,dx)./f5(x);
        pgg_S2 = p1.*integral(fw2, 0, Inf, 'arrayvalued', true);
        
        pg_S2(:,ph)=pgg_S2.';
        ph=ph+1;
        end
    end
end

pot00=reshape(pot0,1,[]);
pot00=pot00.*ones(length(freq),ncondut^2);

       
pot_tot=pot00+pg_S2*1e-3;


pot_tot_matrix=reshape(pot_tot.',ncondut,ncondut,[]);

for nn=1:length(freq)
    pot_total=pot_tot_matrix(:,:,nn);
    pinv=inv(pot_total);
    Y_freq_vetor=pinv*(2*pi*epslonzero)*2*pi*freq(nn)*(-1)^0.5;
    
    Y_freq(nn,:)=reshape(Y_freq_vetor,1,[]);
end


for ii = 1:length(freq)
    Ytotal=reshape(Y_freq(ii,:),ncondut,ncondut);
    Y(:,:,ii)=Ytotal;
end

end

