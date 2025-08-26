close all;
delete(allchild(0));
clear all;

% Toolbox library functions
currmfile = mfilename('fullpath');
currPath = currmfile(1:end-length(mfilename()));
WORKDIR='/home/amauri/Documents/ResearchProjects/LineCableLab/';
addpath(WORKDIR);
addpath([WORKDIR 'ZY_OHTL_pul_funs']);
addpath([WORKDIR 'multilayersoil_funs']);
addpath([WORKDIR 'mode_decomp_funs']);
addpath([WORKDIR 'FD_soil_models_funs']);
addpath([WORKDIR 'line_export_funs']);
addpath([WORKDIR 'Laplace_sim_model_update']);
addpath(fullfile(WORKDIR,'line_export_funs','functions'));
addpath(fullfile(WORKDIR,'line_export_funs','vfit3'));

% Model definition
f=50;
FD_flag = 0;
[line_length,ord,soil,h,d,Geom,Nph,Ncab,ph_order,ph_idx]=LineCableData_fun();
% [line_length,ord,soil,h,d,Geom,Nph]=LineData_fun();
soilFD = soilFD_fun(soil,FD_flag,f);

% Zopts.Papadopoulos=true;
%   allZ_pul = Zcalc_fun(f,ord,Nph,soilFD,h,d,Geom,currPath,'Test',Zopts);

tic
P=cell(Ncab);
w=2*pi*f;
sigma_g=soilFD.sigma_g_total(end);
m_g=soilFD.m_g(end);
e_g=soilFD.e_g_total(end);

% For each unique cable index, find the maximum of the concatenated radii (cols 6 and 9)
outermost_radii = arrayfun(@(cable_num) ...
    max([Geom(Geom(:,1) == cable_num, 6), Geom(Geom(:,1) == cable_num, 9)], [], 'all'), ...
    1:Ncab);

% Earth return impedance matrices considering cables only
e0=8.8541878128e-12;
P0=1/(2*pi*e0);
Pg_self = P0.*P_pet_slf_perf(abs(h(ph_idx)),outermost_radii,Ncab);
Pg_mutual = P0.*P_pet_mut_perf(abs(h(ph_idx)),d(ph_idx,ph_idx),Ncab);
% Pollaczek formulation
% Pg_self=P_papad_slf(h(ph_idx),outermost_radii,0*e_g,m_g,sigma_g,f,Ncab,0);
% Pg_mutual=P_papad_mut(h(ph_idx),d(ph_idx,ph_idx),0*e_g,m_g,sigma_g,f,Ncab,0);
Pext = Pg_self+Pg_mutual;
% Pext = zeros(Ncab);

for i=1:Ncab
    cabledata=Geom(Geom(:,1)==i,2:end);
    ncond_cable=size(cabledata,1); %number of conductor layers in cable i
    for j=1:Ncab
        % loop-based impedance matrix for cable i
        Pcable=zeros(ncond_cable);
        if i==j
            % self-
            for k=1:ncond_cable
                for l=k:ncond_cable
                    Pkl = 0;
                    for m=l:ncond_cable
                        if ~any(isnan(cabledata(k,8:end)))
                            radius_in = cabledata(m, 5);    % Inner radius of insulation layer
                            radius_ex = cabledata(m, 8);    % Outer radius of insulation layer
                            eps_r = cabledata(m, 10);       % Relative permittivity of insulation layer

                            % Call the Potcoef_tub function and accumulate along rows
                            Pkl = Pkl + Potcoef_tub(radius_ex, radius_in, eps_r, w);
                        end
                    end
                    Pcable(k,l) = Pkl;
                end
            end
            Pcable = Pcable + triu(Pcable, 1).';  % mirror the upper triangle onto the lower triangle
        end
        P{i,j} = Pcable + ones(ncond_cable).*Pext(i,j);
    end
end

Pfull=cell2mat(P);
Preduced=bundleReduction(ph_order,Pfull);
Yproposed=1i*w*inv(Preduced)

Yatp=[complex(0,7.42516E-08) complex(0,0); ...
      complex(0,0) complex(0,7.42516E-08)]

calculatePercError(Yproposed,Yatp)


toc
% %% Test with Ametani's formulation
% 
% % Define readable variables to make this bruteforcing less painful
% ii=1;
% cabledata=Geom(Geom(:,1)==ii,2:end);
% jj=1;
% r1=cabledata(jj,4); %core inner radius
% r2=cabledata(jj,5); %core outer radius
% sig_core=1/cabledata(jj,6);
% mu=1;
% eps_ins_core=cabledata(jj,10);
% r3=cabledata(jj+1,4); %sheath inner radius = core outer insulation radius
% r4=cabledata(jj+1,5); % sheath outer radius
% sig_sheath=1/cabledata(jj+1,6);
% eps_ins_sheath=cabledata(jj+1,10);
% r5=cabledata(jj+2,4); %armor inner radius = sheath outer insulation radius
% r6=cabledata(jj+2,5); %armor outer radius
% sig_armor=1/cabledata(jj+2,6);
% eps_ins_armor=cabledata(jj+2,10);
% r7=cabledata(jj+2,8); %armor outer insulation radius
% 
% % z_{cs}=z_{11}+z_{12}+z_{2i}
% z11=Zskeffct_outer(r2,r1,sig_core,mu,w) %core outer impedance
% z12=Zinsu_outer(r3,r2,mu,w) %core insulation outer impedance
% z2i=Zskeffct_inner(r4,r3,sig_sheath,mu,w) %sheath inner core impedance
% zcs=z11+z12+z2i
% 
% % z_{sa}=z_{20}+z_{23}+z_{3i}
% z20=Zskeffct_outer(r4,r3,sig_sheath,mu,w) %sheath outer impedance
% z23=Zinsu_outer(r5,r4,mu,w) %sheath insulation outer impedance
% z3i=Zskeffct_inner(r6,r5,sig_armor,mu,w) %armor inner core impedance
% zsa=z20+z23+z3i
% 
% % z_{a4}=z_{30}+z_{34}
% z30=Zskeffct_outer(r6,r5,sig_armor,mu,w) %armor outer impedance
% z34=Zinsu_outer(r7,r6,mu,w) %armor insulation outer impedance
% za4=z30+z34
% 
% % Z_{ccj}=z_{cs}+z_{sa}+z_{a4}-2z_{2m}-2z_{3m}
% z2m=Zskeffct_mutual(r4,r3,sig_sheath,mu,w) %sheath mutual
% z3m=Zskeffct_mutual(r6,r5,sig_armor,mu,w) %armor mutual
% 
% Zcc=zcs+zsa+za4-2*z2m-2*z3m
% Zss=zsa+za4-2*z3m
% Zaa=za4
% Zcs=zsa+za4-z2m-2*z3m
% Zca=za4-z3m
% Zsa=Zca
% 
% Zi=[Zcc Zcs Zca; ...
%     Zcs Zss Zsa; ...
%     Zca Zsa Zaa]

% Z =
% 
% [Z11, Z12,   0]
% [Z21, Z22, Z23]
% [  0, Z32, Z33]
% 
% Transformed symbolic impedance matrix Z_prime:
% [Z11 + Z12 + Z21 + Z22 + Z23 + Z32 + Z33, Z12 + Z22 + Z23 + Z32 + Z33, Z23 + Z33]
% [            Z21 + Z22 + Z23 + Z32 + Z33,       Z22 + Z23 + Z32 + Z33, Z23 + Z33]
% [                              Z32 + Z33,                   Z32 + Z33,       Z33]
% 
