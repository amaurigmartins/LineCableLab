close all
% clear all
clc

currmfile = mfilename('fullpath');
currPath = currmfile(1:end-length(mfilename()));
addpath([currPath 'ZY_OHTL_pul_funs']);
addpath([currPath 'mode_decomp_funs']);
addpath([currPath 'FD_soil_models_funs']);

points_per_dec=100;
f_dec1=1:10/points_per_dec:10;
f_dec2=10:100/points_per_dec:100;
f_dec3=100:1000/(points_per_dec):1000;
f_dec4=1000:9000/(points_per_dec):9000; %%%%%%%%%%%%%%%%%% this decade
f_dec5=9000:100000/(points_per_dec/10):100000;
f_dec6=100000:1000000/points_per_dec:1000000;
f=transpose([f_dec1(1:length(f_dec1)-1) f_dec2(1:length(f_dec2)-1) f_dec3(1:length(f_dec3)-1) f_dec4(1:length(f_dec4)-1) f_dec5(1:length(f_dec5)-1) f_dec6]);

num_modes = 6;

const_model_dataset={ ...
    'FD0_Ztot_Carson_Ytot_Imag_rho100_eps15';...
    'FD0_Ztot_Carson_Ytot_Imag_rho1000_eps5';...
    'FD0_Ztot_Carson_Ytot_Imag_rho5000_eps3';...
    'FD0_Ztot_Wise_Ytot_Wise_rho100_eps15';...
    'FD0_Ztot_Wise_Ytot_Wise_rho1000_eps5';...
    'FD0_Ztot_Wise_Ytot_Wise_rho5000_eps3'};

LS_model_dataset={ ...
    'FD1_Ztot_Carson_Ytot_Imag_rho104_eps15';...
    'FD1_Ztot_Carson_Ytot_Imag_rho1097_eps5';...
    'FD1_Ztot_Carson_Ytot_Imag_rho5790_eps3';...
    'FD1_Ztot_Wise_Ytot_Wise_rho104_eps15';...
    'FD1_Ztot_Wise_Ytot_Wise_rho1097_eps5';...
    'FD1_Ztot_Wise_Ytot_Wise_rho5790_eps3'};

cigre_model_dataset={ ...
    'FD9_Ztot_Carson_Ytot_Imag_rho100_eps15';...
    'FD9_Ztot_Carson_Ytot_Imag_rho1000_eps5';...
    'FD9_Ztot_Carson_Ytot_Imag_rho5000_eps3';...
    'FD9_Ztot_Wise_Ytot_Wise_rho100_eps15';...
    'FD9_Ztot_Wise_Ytot_Wise_rho1000_eps5';...
    'FD9_Ztot_Wise_Ytot_Wise_rho5000_eps3'};


%% plots for a specific mode
mode = 3;
ds_name='const_model';
this_dataset=eval([ds_name '_dataset']);
rho_names={'100', '1000', '5000'};


set_plot_params()
figname=['WisCarPropag_noratio___' ds_name '_mode' num2str(mode)]
figure('Name',figname)

subplot(2,1,1)
for idx=1:3
a_wise=read_from_file(this_dataset{idx+3},'a_dis');
semilogx(f,a_wise(:,mode),'DisplayName', ['Wise - \rho = ' rho_names{idx} ' \Omega.m'], 'LineStyle','-');hold all

a_carson=read_from_file(this_dataset{idx},'a_dis');
semilogx(f,a_carson(:,mode),'DisplayName', ['Carson - \rho = ' rho_names{idx} ' \Omega.m'], 'LineStyle',':');hold all

end
title(['Attenuation constant - mode #' num2str(mode) ' (a)'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha [Np/m]')
grid on
legend

subplot(2,1,2)
for idx=1:3
vel_wise=read_from_file(this_dataset{idx+3},'vel_dis');
semilogx(f,vel_wise(:,mode),'DisplayName', ['Wise - \rho = ' rho_names{idx} ' \Omega.m'], 'LineStyle','-');hold all
vel_carson=read_from_file(this_dataset{idx},'vel_dis');
semilogx(f,vel_carson(:,mode),'DisplayName', ['Carson - \rho = ' rho_names{idx} ' \Omega.m'], 'LineStyle',':');hold all
end
title(['Phase velocity - mode #' num2str(mode) ' (b)'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon [m/s]')
grid on
legend


%% plots for all modes

idx=3; % 1 = 100, 2 = 1000, 3 = 5000
a_carson=read_from_file(this_dataset{idx},'a_dis');
vel_carson=read_from_file(this_dataset{idx},'vel_dis');
a_wise=read_from_file(this_dataset{idx+3},'a_dis');
vel_wise=read_from_file(this_dataset{idx+3},'vel_dis');

set_plot_params([1 2 2 1 1 2 1])
figure('Name',['WisCar_AtnConst_noratio___' ds_name '_rho' rho_names{idx}])

for i=1:num_modes
    subplot(2,3,i)
    semilogx(f,a_carson(:,i),'DisplayName', ['Carson']);hold all
    semilogx(f,a_wise(:,i),'DisplayName', ['Wise']);hold all
    title(['mode #' num2str(i)])
    axis tight
    xlabel('Frequency [Hz]')
    ylabel('\alpha [Np/m]')
    grid on
    legend
end

set_plot_params([1 2 2 1 1 2 1])
figure('Name',['WisCar_PhaseVel_noratio_' ds_name '_' rho_names{idx}])

for i=1:num_modes
    subplot(2,3,i)
    semilogx(f,vel_carson(:,i),'DisplayName', ['Carson']);hold all
    semilogx(f,vel_wise(:,i),'DisplayName', ['Wise']);hold all
    title(['mode #' num2str(i)])
    axis tight
    xlabel('Frequency [Hz]')
    ylabel('\upsilon [m/s]')
    grid on
    legend
end
