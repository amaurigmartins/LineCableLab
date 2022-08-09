close all
clear all
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


Ztot_Carson_Ytot_Imag_const_model_dataset={ ...
    'FD0_Ztot_Carson_Ytot_Imag_rho100_eps15';...
    'FD0_Ztot_Carson_Ytot_Imag_rho1000_eps5';...
    'FD0_Ztot_Carson_Ytot_Imag_rho5000_eps3'};

Ztot_Carson_Ytot_Imag_LS_model_dataset={ ...
    'FD1_Ztot_Carson_Ytot_Imag_rho104_eps15';...
    'FD1_Ztot_Carson_Ytot_Imag_rho1097_eps5';...
    'FD1_Ztot_Carson_Ytot_Imag_rho5790_eps3'};

Ztot_Carson_Ytot_Imag_cigre_model_dataset={ ...
    'FD9_Ztot_Carson_Ytot_Imag_rho100_eps15';...
    'FD9_Ztot_Carson_Ytot_Imag_rho1000_eps5';...
    'FD9_Ztot_Carson_Ytot_Imag_rho5000_eps3'};




Ztot_Wise_Ytot_Wise_const_model_dataset={ ...
    'FD0_Ztot_Wise_Ytot_Wise_rho100_eps15';...
    'FD0_Ztot_Wise_Ytot_Wise_rho1000_eps5';...
    'FD0_Ztot_Wise_Ytot_Wise_rho5000_eps3'};

Ztot_Wise_Ytot_Wise_LS_model_dataset={ ...
    'FD1_Ztot_Wise_Ytot_Wise_rho104_eps15';...
    'FD1_Ztot_Wise_Ytot_Wise_rho1097_eps5';...
    'FD1_Ztot_Wise_Ytot_Wise_rho5790_eps3'};

Ztot_Wise_Ytot_Wise_cigre_model_dataset={ ...
    'FD9_Ztot_Wise_Ytot_Wise_rho100_eps15';...
    'FD9_Ztot_Wise_Ytot_Wise_rho1000_eps5';...
    'FD9_Ztot_Wise_Ytot_Wise_rho5000_eps3'};

ds_name={'const_model','LS_model', 'cigre_model'};
% this_dataset=eval([ds_name '_dataset']);
idx=3; % 1 = 100, 2 = 1000, 3 = 10000
rho_names={'100', '1000', '5000'};
formula={'Ztot_Carson_Ytot_Imag', 'Ztot_Wise_Ytot_Wise'};

set_plot_params([1 2 2 1 1 2 1])
figure('Name',['FD_CP_CharImpMag_noratio___'])

for i=1:num_modes
    subplot(2,3,i)
    for j=1:2
        for k=1:3
            zch=read_from_file(eval([formula{j} '_' ds_name{k} '_dataset{idx}']),'Zch_mod');
            semilogx(f,abs(zch(:,i)),'DisplayName', [formula{j}(6:9) '__' ds_name{k}(1:end-6) '__rho' rho_names{k}]);hold all
        end
    end
    title(['mode #' num2str(i)])
    axis tight
    xlabel('Frequency [Hz]')
    ylabel('Z_{ch} magnitude [\Omega/m]')
    grid on
    legend
end

set_plot_params([1 2 2 1 1 2 1])
figure('Name',['FD_CP_CharImpMag_noratio___'])

for i=1:num_modes
    subplot(2,3,i)
    for j=1:2
        for k=1:3
            zch=read_from_file(eval([formula{j} '_' ds_name{k} '_dataset{idx}']),'Zch_mod');
            semilogx(f,radtodeg(angle(zch(:,i))),'DisplayName', [formula{j}(6:9) '__' ds_name{k}(1:end-6) '__rho' rho_names{k}]);hold all
        end
    end
    title(['mode #' num2str(i)])
    axis tight
    xlabel('Frequency [Hz]')
    ylabel('Z_{ch} angle [deg]')
    grid on
    legend
end





set_plot_params([1 2 2 1 1 2 1])
figure('Name',['FD_CP_AtnConst_noratio___'])

for i=1:num_modes
    subplot(2,3,i)
    for j=1:2
        for k=1:3
            a=read_from_file(eval([formula{j} '_' ds_name{k} '_dataset{idx}']),'a_dis');
            semilogx(f,a(:,i),'DisplayName', [formula{j}(6:9) '__' ds_name{k}(1:end-6) '__rho' rho_names{k}]);hold all
        end
    end
    title(['mode #' num2str(i)])
    axis tight
    xlabel('Frequency [Hz]')
    ylabel('\alpha [Np/m]')
    grid on
    legend
end








set_plot_params([1 2 2 1 1 2 1])
figure('Name',['FD_CP_PhaseVel_noratio___'])

for i=1:num_modes
    subplot(2,3,i)
    for j=1:2
        for k=1:3
            vel=read_from_file(eval([formula{j} '_' ds_name{k} '_dataset{idx}']),'vel_dis');
            semilogx(f,vel(:,i),'DisplayName', [formula{j}(6:9) '__' ds_name{k}(1:end-6) '__rho' rho_names{k}]);hold all
        end
    end
    title(['mode #' num2str(i)])
    axis tight
    xlabel('Frequency [Hz]')
    ylabel('\upsilon [m/s]')
    grid on
    legend
end




