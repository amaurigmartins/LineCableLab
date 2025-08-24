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
Z=cell(Ncab);
w=2*pi*f;
sigma_g=soilFD.sigma_g_total(end);
m_g=soilFD.m_g(end);
e_g=soilFD.e_g_total(end);

% For each unique cable index, find the maximum of the concatenated radii (cols 6 and 9)
outermost_radii = arrayfun(@(cable_num) ...
    max([Geom(Geom(:,1) == cable_num, 6), Geom(Geom(:,1) == cable_num, 9)], [], 'all'), ...
    1:Ncab);

% Earth return impedance matrices considering cables only 
Zg_self=Z_papad_slf(h(ph_idx),outermost_radii,0*e_g ,m_g,sigma_g,f,Ncab,0);
Zg_mutual=Z_papad_mut(h(ph_idx),d(ph_idx,ph_idx),0*e_g,m_g,sigma_g,f,Ncab,0);
Zext = Zg_self+Zg_mutual;
% Zext=zeros(Ncab)

for i=1:Ncab
    cabledata=Geom(Geom(:,1)==i,2:end);
    ncond_cable=size(cabledata,1); %number of conductor layers in cable i
    for j=1:Ncab
        % loop-based impedance matrix for cable i
        Zloop=zeros(ncond_cable);
        if i==j
            % self-
            for k=1:ncond_cable % banded matrix
                % this layer outer surface impedance 
                rin_thislayer=cabledata(k,4);
                rext_thislayer=cabledata(k,5);
                sig_c_thislayer=1/cabledata(k,6);
                mu_c_thislayer=cabledata(k,7);
                Zouter_thislayer=Zskeffct_outer(rext_thislayer,rin_thislayer,sig_c_thislayer,mu_c_thislayer,w);
                % this layer insulation impedance (if present)
                if ~any(isnan(cabledata(k,8:end)))
                    rext_thisinsu=cabledata(k,8);
                    mu_thisinsu=cabledata(k,9);
                    Zinsu_thislayer=Zinsu_outer(rext_thisinsu,rext_thislayer,mu_thisinsu,w);
                else
                    Zinsu_thislayer=0;
                end
                % next layer inner surface impedance
                if k<ncond_cable
                    rin_nextlayer=cabledata(k+1,4);
                    rext_nextlayer=cabledata(k+1,5);
                    sig_c_nextlayer=1/cabledata(k+1,6);
                    mu_c_nextlayer=cabledata(k+1,7);
                    Zinner_nextlayer=Zskeffct_inner(rext_nextlayer,rin_nextlayer,sig_c_nextlayer,mu_c_nextlayer,w);
                    Zmutual=Zskeffct_mutual(rext_nextlayer,rin_nextlayer,sig_c_nextlayer,mu_c_nextlayer,w);
                else
                    Zinner_nextlayer=Zext(i,j); % self- earth return impedance (underground conductor)
                    Zmutual=0;
                end
                Zloop(k,k)=Zouter_thislayer+Zinsu_thislayer+Zinner_nextlayer;
                if k<ncond_cable
                    Zloop(k,k+1)=-Zmutual;
                    Zloop(k+1,k)=-Zmutual;
                end
            end
        else
            % mutual-
            Zloop(end,end)=Zext(i,j);
        end
        Z{i,j}=transformZloop(Zloop);
    end
end

Zfull=cell2mat(Z);
Zproposed=bundleReduction(ph_order,Zfull)
Zatp=[complex(1.01276E-04,6.07996E-05) complex(5.51412E-07,-4.06385E-06); ...
      complex(5.51412E-07,-4.06385E-06) complex(1.01276E-04,6.07996E-05)]

calculatePercError(Zproposed,Zatp)


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
