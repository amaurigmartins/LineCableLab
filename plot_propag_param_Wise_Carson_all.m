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

%% Wise x Carson

% **************************** 100 ohm.m, contant soil model
a_ref=read_from_file(const_model_dataset{1},'a_dis');
vel_ref=read_from_file(const_model_dataset{1},'vel_dis');
a=read_from_file(const_model_dataset{4},'a_dis');
vel=read_from_file(const_model_dataset{4},'vel_dis');
yya = a ./ a_ref;
yyv = vel ./ vel_ref;

set_plot_params([1 2 2 1 1 2 1])
figure('Name','WiseCarson___rho___const')
subplot(2,3,1)
for i=1:num_modes
    semilogx(f,yya(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['\rho = 100 \Omega.m, ' char(949) '_{r} = 15'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(2,3,4)
for i=1:num_modes
    semilogx(f,yyv(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['\rho = 100 \Omega.m, ' char(949) '_{r} = 15'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

% 1000 ohm.m, contant soil model
a_ref=read_from_file(const_model_dataset{2},'a_dis');
vel_ref=read_from_file(const_model_dataset{2},'vel_dis');
a=read_from_file(const_model_dataset{5},'a_dis');
vel=read_from_file(const_model_dataset{5},'vel_dis');
yya = a ./ a_ref;
yyv = vel ./ vel_ref;

subplot(2,3,2)
for i=1:num_modes
    semilogx(f,yya(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['\rho = 1000 \Omega.m, ' char(949) '_{r} = 5'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(2,3,5)
for i=1:num_modes
    semilogx(f,yyv(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['\rho = 1000 \Omega.m, ' char(949) '_{r} = 5'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

% 10000 ohm.m, contant soil model
a_ref=read_from_file(const_model_dataset{3},'a_dis');
vel_ref=read_from_file(const_model_dataset{3},'vel_dis');
a=read_from_file(const_model_dataset{6},'a_dis');
vel=read_from_file(const_model_dataset{6},'vel_dis');
yya = a ./ a_ref;
yyv = vel ./ vel_ref;

subplot(2,3,3)
for i=1:num_modes
    semilogx(f,yya(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['\rho = 5000 \Omega.m, ' char(949) '_{r} = 3'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(2,3,6)
for i=1:num_modes
    semilogx(f,yyv(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['\rho = 5000 \Omega.m, ' char(949) '_{r} = 3'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

% **************************** 100 ohm.m, LS soil model
a_ref=read_from_file(LS_model_dataset{1},'a_dis');
vel_ref=read_from_file(LS_model_dataset{1},'vel_dis');
a=read_from_file(LS_model_dataset{4},'a_dis');
vel=read_from_file(LS_model_dataset{4},'vel_dis');
yya = a ./ a_ref;
yyv = vel ./ vel_ref;

set_plot_params([1 2 2 1 1 2 1])
figure('Name','WiseCarson___rho___LS')
subplot(2,3,1)
for i=1:num_modes
    semilogx(f,yya(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['\rho = 100 \Omega.m, ' char(949) '_{r} = 15'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(2,3,4)
for i=1:num_modes
    semilogx(f,yyv(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['\rho = 100 \Omega.m, ' char(949) '_{r} = 15'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

% 1000 ohm.m, LS soil model
a_ref=read_from_file(LS_model_dataset{2},'a_dis');
vel_ref=read_from_file(LS_model_dataset{2},'vel_dis');
a=read_from_file(LS_model_dataset{5},'a_dis');
vel=read_from_file(LS_model_dataset{5},'vel_dis');
yya = a ./ a_ref;
yyv = vel ./ vel_ref;

subplot(2,3,2)
for i=1:num_modes
    semilogx(f,yya(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['\rho = 1000 \Omega.m, ' char(949) '_{r} = 5'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(2,3,5)
for i=1:num_modes
    semilogx(f,yyv(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['\rho = 1000 \Omega.m, ' char(949) '_{r} = 5'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

% 10000 ohm.m, LS soil model
a_ref=read_from_file(LS_model_dataset{3},'a_dis');
vel_ref=read_from_file(LS_model_dataset{3},'vel_dis');
a=read_from_file(LS_model_dataset{6},'a_dis');
vel=read_from_file(LS_model_dataset{6},'vel_dis');
yya = a ./ a_ref;
yyv = vel ./ vel_ref;

subplot(2,3,3)
for i=1:num_modes
    semilogx(f,yya(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['\rho = 5000 \Omega.m, ' char(949) '_{r} = 3'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(2,3,6)
for i=1:num_modes
    semilogx(f,yyv(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['\rho = 5000 \Omega.m, ' char(949) '_{r} = 3'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

% **************************** 100 ohm.m, cigre soil model
a_ref=read_from_file(cigre_model_dataset{1},'a_dis');
vel_ref=read_from_file(cigre_model_dataset{1},'vel_dis');
a=read_from_file(cigre_model_dataset{4},'a_dis');
vel=read_from_file(cigre_model_dataset{4},'vel_dis');
yya = a ./ a_ref;
yyv = vel ./ vel_ref;

set_plot_params([1 2 2 1 1 2 1])
figure('Name','WiseCarson___rho___cigre')
subplot(2,3,1)
for i=1:num_modes
    semilogx(f,yya(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['\rho = 100 \Omega.m, ' char(949) '_{r} = 15'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(2,3,4)
for i=1:num_modes
    semilogx(f,yyv(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['\rho = 100 \Omega.m, ' char(949) '_{r} = 15'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

% 1000 ohm.m, cigre soil model
a_ref=read_from_file(cigre_model_dataset{2},'a_dis');
vel_ref=read_from_file(cigre_model_dataset{2},'vel_dis');
a=read_from_file(cigre_model_dataset{5},'a_dis');
vel=read_from_file(cigre_model_dataset{5},'vel_dis');
yya = a ./ a_ref;
yyv = vel ./ vel_ref;

subplot(2,3,2)
for i=1:num_modes
    semilogx(f,yya(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['\rho = 1000 \Omega.m, ' char(949) '_{r} = 5'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(2,3,5)
for i=1:num_modes
    semilogx(f,yyv(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['\rho = 1000 \Omega.m, ' char(949) '_{r} = 5'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

% 10000 ohm.m, cigre soil model
a_ref=read_from_file(cigre_model_dataset{3},'a_dis');
vel_ref=read_from_file(cigre_model_dataset{3},'vel_dis');
a=read_from_file(cigre_model_dataset{6},'a_dis');
vel=read_from_file(cigre_model_dataset{6},'vel_dis');
yya = a ./ a_ref;
yyv = vel ./ vel_ref;

subplot(2,3,3)
for i=1:num_modes
    semilogx(f,yya(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['\rho = 5000 \Omega.m, ' char(949) '_{r} = 3'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(2,3,6)
for i=1:num_modes
    semilogx(f,yyv(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['\rho = 5000 \Omega.m, ' char(949) '_{r} = 3'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

% **************************** 100 ohm.m, const soil model
a_ref=read_from_file(const_model_dataset{1},'a_dis');
vel_ref=read_from_file(const_model_dataset{1},'vel_dis');
a=read_from_file(const_model_dataset{4},'a_dis');
vel=read_from_file(const_model_dataset{4},'vel_dis');
yya = a ./ a_ref;
yyv = vel ./ vel_ref;

set_plot_params([1 2 2 1 1 2 1])
figure('Name','WiseCarson___FDmodel___rho100')
subplot(2,3,1)
for i=1:num_modes
    semilogx(f,yya(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['Constant soil parameters, \rho = 100 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(2,3,4)
for i=1:num_modes
    semilogx(f,yyv(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['Constant soil parameters, \rho = 100 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

% 100 ohm.m, LS soil model
a_ref=read_from_file(LS_model_dataset{1},'a_dis');
vel_ref=read_from_file(LS_model_dataset{1},'vel_dis');
a=read_from_file(LS_model_dataset{4},'a_dis');
vel=read_from_file(LS_model_dataset{4},'vel_dis');
yya = a ./ a_ref;
yyv = vel ./ vel_ref;

subplot(2,3,2)
for i=1:num_modes
    semilogx(f,yya(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['LS soil model, \rho = 100 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(2,3,5)
for i=1:num_modes
    semilogx(f,yyv(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['LS soil model, \rho = 100 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

% 100 ohm.m, Cigre  soil model
a_ref=read_from_file(cigre_model_dataset{1},'a_dis');
vel_ref=read_from_file(cigre_model_dataset{1},'vel_dis');
a=read_from_file(cigre_model_dataset{4},'a_dis');
vel=read_from_file(cigre_model_dataset{4},'vel_dis');
yya = a ./ a_ref;
yyv = vel ./ vel_ref;

subplot(2,3,3)
for i=1:num_modes
    semilogx(f,yya(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['Cigre soil model, \rho = 100 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(2,3,6)
for i=1:num_modes
    semilogx(f,yyv(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['Cigre soil model, \rho = 100 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

% **************************** 1000 ohm.m, const soil model
a_ref=read_from_file(const_model_dataset{2},'a_dis');
vel_ref=read_from_file(const_model_dataset{2},'vel_dis');
a=read_from_file(const_model_dataset{5},'a_dis');
vel=read_from_file(const_model_dataset{5},'vel_dis');
yya = a ./ a_ref;
yyv = vel ./ vel_ref;

set_plot_params([1 2 2 1 1 2 1])
figure('Name','WiseCarson___FDmodel___rho1000')
subplot(2,3,1)
for i=1:num_modes
    semilogx(f,yya(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['Constant soil parameters, \rho = 1000 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(2,3,4)
for i=1:num_modes
    semilogx(f,yyv(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['Constant soil parameters, \rho = 1000 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

% 1000 ohm.m, LS soil model
a_ref=read_from_file(LS_model_dataset{2},'a_dis');
vel_ref=read_from_file(LS_model_dataset{2},'vel_dis');
a=read_from_file(LS_model_dataset{5},'a_dis');
vel=read_from_file(LS_model_dataset{5},'vel_dis');
yya = a ./ a_ref;
yyv = vel ./ vel_ref;

subplot(2,3,2)
for i=1:num_modes
    semilogx(f,yya(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['LS soil model, \rho = 1000 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(2,3,5)
for i=1:num_modes
    semilogx(f,yyv(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['LS soil model, \rho = 1000 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

% 1000 ohm.m, Cigre  soil model
a_ref=read_from_file(cigre_model_dataset{2},'a_dis');
vel_ref=read_from_file(cigre_model_dataset{2},'vel_dis');
a=read_from_file(cigre_model_dataset{5},'a_dis');
vel=read_from_file(cigre_model_dataset{5},'vel_dis');
yya = a ./ a_ref;
yyv = vel ./ vel_ref;

subplot(2,3,3)
for i=1:num_modes
    semilogx(f,yya(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['Cigre soil model, \rho = 1000 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(2,3,6)
for i=1:num_modes
    semilogx(f,yyv(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['Cigre soil model, \rho = 1000 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

% **************************** 10000 ohm.m, const soil model
a_ref=read_from_file(const_model_dataset{3},'a_dis');
vel_ref=read_from_file(const_model_dataset{3},'vel_dis');
a=read_from_file(const_model_dataset{6},'a_dis');
vel=read_from_file(const_model_dataset{6},'vel_dis');
yya = a ./ a_ref;
yyv = vel ./ vel_ref;

set_plot_params([1 2 2 1 1 2 1])
figure('Name','WiseCarson___FDmodel___rho10000')
subplot(2,3,1)
for i=1:num_modes
    semilogx(f,yya(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['Constant soil parameters, \rho = 5000 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(2,3,4)
for i=1:num_modes
    semilogx(f,yyv(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['Constant soil parameters, \rho = 5000 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

% 10000 ohm.m, LS soil model
a_ref=read_from_file(LS_model_dataset{3},'a_dis');
vel_ref=read_from_file(LS_model_dataset{3},'vel_dis');
a=read_from_file(LS_model_dataset{6},'a_dis');
vel=read_from_file(LS_model_dataset{6},'vel_dis');
yya = a ./ a_ref;
yyv = vel ./ vel_ref;

subplot(2,3,2)
for i=1:num_modes
    semilogx(f,yya(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['LS soil model, \rho = 5000 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(2,3,5)
for i=1:num_modes
    semilogx(f,yyv(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['LS soil model, \rho = 5000 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

% 10000 ohm.m, Cigre  soil model
a_ref=read_from_file(cigre_model_dataset{3},'a_dis');
vel_ref=read_from_file(cigre_model_dataset{3},'vel_dis');
a=read_from_file(cigre_model_dataset{6},'a_dis');
vel=read_from_file(cigre_model_dataset{6},'vel_dis');
yya = a ./ a_ref;
yyv = vel ./ vel_ref;

subplot(2,3,3)
for i=1:num_modes
    semilogx(f,yya(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['Cigre soil model, \rho = 5000 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(2,3,6)
for i=1:num_modes
    semilogx(f,yyv(:,i),'DisplayName', ['mode #' num2str(i)]);hold all
end
title(['Cigre soil model, \rho = 5000 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

