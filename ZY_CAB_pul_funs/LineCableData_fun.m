function [length,Ncon,soil,h,d,Geom,Nph,Ncables,ph_order,ph_idx]=LineCableData_fun()
% Line geometry
% 1 column  -- number of phase (set to 0 for Kron reduction)
% 2 column  -- x position of conduntor center in meters
% 3 column  -- y position of coductor center in meters
% 4 column  -- internal radius of conductor
% 5 column  -- external radius of conductor
% 6 column  -- resistivity of conductor
% 7 column  -- permeability of conductor
% 8 column  -- external radius of conductor insulation
% 9 column  -- relative permeability of conductor insulation
% 10 column -- relative permittivity of conductor insulation
% 11 column -- external radius of conductor insulation

length = 1000.000000;

% cross section data from Papagiannis et al
% conductor properties from Gustavsen et al
h=-0.75;
xa=-0.15;
xb=0;
xc=0.15;
mu_r=1;
r1=0;
r2=0.0127; %core
rho_c=3.4643e-8; 
r3=0.0228; %insulation 1
er1=2.486; 
r4=0.0254; %sheath
rho_s=1.718e-8; 
r5=0.0279; %insulation 2
er2=2.856; 
r6=r5+0.5*(r4-r3); % armor w/ foil thickness assumed half of the sheath
rho_a=2.8264e-8; 
r7=r6+(r5-r4); %insulation 3
er3=2.4;

% --- alternative, maybe better
%      circuit    | geom    | conductor                 | insulation
%    1----2-------3----4----5----6-----7--------8-------9-----10-----11 
%   cabID phID    horz vert r_in r_ext rho_cond mu_cond r_ins mu_ins eps_ins
Geom=[ ...
     1    1       xa   h    r1   r2    rho_c    mu_r    r3    mu_r   er1
     1    0       nan  nan  r3   r4    rho_s    mu_r    r5    mu_r   er2
     1    0       nan  nan  r5   r6    rho_a    mu_r    r7    mu_r   er3
     2    2       xb   h    r1   r2    rho_c    mu_r    r3    mu_r   er1
     2    0       nan  nan  r3   r4    rho_s    mu_r    r5    mu_r   er2
     2    0       nan  nan  r5   r6    rho_a    mu_r    r7    mu_r   er3
     3    3       xc   h    r1   r2    rho_c    mu_r    r3    mu_r   er1
     3    0       nan  nan  r3   r4    rho_s    mu_r    r5    mu_r   er2
     3    0       nan  nan  r5   r6    rho_a    mu_r    r7    mu_r   er3
];

Ncables=size(unique(Geom(:,1)),1);
Ncon = size(Geom,1);

% Variables
e0=8.854187817e-12;
m0=4*pi*1e-7;
% Rearrange data matrix to perform bundle reduction
ph_order=Geom(:,2);
ph_idx=find(~(isnan(Geom(:,3)) | isnan(Geom(:,4))));
% Geom = reorderGeoMatrix(Geom);
% ph_reorder=Geom(:,2);
Nph = unique(ph_order);
Nph = max(Nph(Nph~=0));
% Height of line calculation
[h]=height_fun(Ncon,Geom);
% Distances calculations
[d]=distance_fun(Ncon,Geom);
% Earth electric parameters
soil.erg = 15.000000;
soil.m_g = m0*1.000000;
soil.sigma_g = 1/100.000000;
end

%% reorderGeoMatrix
function MF = reorderGeoMatrix(M)
M = sortrows(M,2);

zero_ph_row = find(M(:,2) == 0);
if ~isempty(zero_ph_row)
     M0 = M(zero_ph_row(1):zero_ph_row(end),:);
     M1 = M(zero_ph_row(end)+1:end,:);
     PH = unique(M(zero_ph_row(end)+1:end,1));
     
else
    M1 = M;
    M0 = [];
end

MF = [];
k = 1;
i = 1;
N = size(M1,2);

Nph = unique(M1(:,2));
Nph = size(Nph(Nph~=0),1);

for i = 1:2
    for j = 1:Nph
        if i == 1
            ph_row = find(M1(:,1)==j);
            if ~isempty(ph_row)
                MF(j,:) = M1(ph_row(1),:);
                M1(ph_row(1),:) = [];
            end
        else
            ph_row = find(M1(:,1)==j);
            if ~isempty(ph_row)
                MF_aux = M1(ph_row,:);
                M1(ph_row,:) = [];
                MF = [MF;MF_aux];
            end
        end
    end
end
           
MF = [MF; M0];

end

%% height_fun
% Set the conductors at the same height 
function [h]=height_fun(num,Geom)
h=zeros(1,num);

for k=1:1:num
    h(1,k)=Geom(k,4);
end

end

%% distance_fun
%Intermediate distance between consecutive conductors
function [d]=distance_fun(num,Geom)

%create a matrix where we will save the distance between the different
%conductors

d=zeros(num,num);

int_cond=1;

for x=1:1:num

    for y=1:1:num
        if x==y
            d(x,y)=abs(Geom(x,3)-Geom(y,3));
            int_cond=1;

        elseif abs(x-y)==int_cond
            d(x,y)=abs(Geom(x,3)-Geom(y,3));
            d(y,x)=d(x,y);
            int_cond=int_cond+1;

        end

    end
    int_cond=1;

end
end