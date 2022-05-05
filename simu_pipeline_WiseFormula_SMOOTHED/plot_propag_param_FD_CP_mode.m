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


% Mode #1 - Blue: Travels through conductor 5 (ground wire) - Returns through conductor 6 (ground wire)
% Mode #2 - Red: Travels through conductors 5 and 6 (ground wires) - Returns through conductors 1, 2, 3 (line phases)
% Mode #3 - Orange: Travels through all metallic conductors - Returns through infinite earth - This is the pure ground mode
% Mode #4 - Purple: Travels through conductors 1, 2, 3, 5, 6 (line phases and ground wires) - Returns through conductor 4 (pipe)
% Mode #5 - Green: Travels through conductor 1 (line phase) - Returns through conductor 3 (line phase)
% Mode #6 - Cyan: Travels through conductor 2 (line phase) - Returns through conductors 1, 3 (line phases)

mode = 3;

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

%% FD x CP for a specific mode

set_plot_params([1 2 2 1 1 4 1])
figure('Name',['FDCP___rho___LS_mode' num2str(mode)])
subplot(1,2,1)

a_ref=read_from_file(const_model_dataset{4},'a_dis');
a=read_from_file(LS_model_dataset{4},'a_dis');
yya = a ./ a_ref;
semilogx(f,yya(:,mode),'DisplayName', '\rho = 100 \Omega.m');hold all

a_ref=read_from_file(const_model_dataset{5},'a_dis');
a=read_from_file(LS_model_dataset{5},'a_dis');
yya = a ./ a_ref;
semilogx(f,yya(:,mode),'DisplayName', '\rho = 1000 \Omega.m');hold all

a_ref=read_from_file(const_model_dataset{6},'a_dis');
a=read_from_file(LS_model_dataset{6},'a_dis');
yya = a ./ a_ref;
semilogx(f,yya(:,mode),'DisplayName', '\rho = 5000 \Omega.m');hold all

title(['mode # ' num2str(mode) ' - Attenuation constant, LS model'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(1,2,2)

vel_ref=read_from_file(const_model_dataset{4},'vel_dis');
vel=read_from_file(LS_model_dataset{4},'vel_dis');
yyv = vel ./ vel_ref;
semilogx(f,yyv(:,mode),'DisplayName', '\rho = 100 \Omega.m');hold all

vel_ref=read_from_file(const_model_dataset{5},'vel_dis');
vel=read_from_file(LS_model_dataset{5},'vel_dis');
yyv = vel ./ vel_ref;
semilogx(f,yyv(:,mode),'DisplayName', '\rho = 1000 \Omega.m');hold all

vel_ref=read_from_file(const_model_dataset{6},'vel_dis');
vel=read_from_file(LS_model_dataset{6},'vel_dis');
yyv = vel ./ vel_ref;
semilogx(f,yyv(:,mode),'DisplayName', '\rho = 5000 \Omega.m');hold all

title(['mode # ' num2str(mode) ' - Phase velocity, LS model'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

set_plot_params([1 2 2 1 1 4 1])
figure('Name',['FDCP___rho___cigre_mode' num2str(mode)])
subplot(1,2,1)

a_ref=read_from_file(const_model_dataset{4},'a_dis');
a=read_from_file(cigre_model_dataset{4},'a_dis');
yya = a ./ a_ref;
semilogx(f,yya(:,mode),'DisplayName', '\rho = 100 \Omega.m');hold all

a_ref=read_from_file(const_model_dataset{5},'a_dis');
a=read_from_file(cigre_model_dataset{5},'a_dis');
yya = a ./ a_ref;
semilogx(f,yya(:,mode),'DisplayName', '\rho = 1000 \Omega.m');hold all

a_ref=read_from_file(const_model_dataset{6},'a_dis');
a=read_from_file(cigre_model_dataset{6},'a_dis');
yya = a ./ a_ref;
semilogx(f,yya(:,mode),'DisplayName', '\rho = 5000 \Omega.m');hold all

title(['mode # ' num2str(mode) ' - Attenuation constant, Cigre model'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(1,2,2)

vel_ref=read_from_file(const_model_dataset{4},'vel_dis');
vel=read_from_file(cigre_model_dataset{4},'vel_dis');
yyv = vel ./ vel_ref;
semilogx(f,yyv(:,mode),'DisplayName', '\rho = 100 \Omega.m');hold all

vel_ref=read_from_file(const_model_dataset{5},'vel_dis');
vel=read_from_file(cigre_model_dataset{5},'vel_dis');
yyv = vel ./ vel_ref;
semilogx(f,yyv(:,mode),'DisplayName', '\rho = 1000 \Omega.m');hold all

vel_ref=read_from_file(const_model_dataset{6},'vel_dis');
vel=read_from_file(cigre_model_dataset{6},'vel_dis');
yyv = vel ./ vel_ref;
semilogx(f,yyv(:,mode),'DisplayName', '\rho = 5000 \Omega.m');hold all

title(['mode # ' num2str(mode) ' - Phase velocity, Cigre model'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend


set_plot_params([1 2 2 1 1 4 1])
figure('Name',['FDCP___FDmodel___rho100_mode' num2str(mode)])
subplot(1,2,1)

a_ref=read_from_file(const_model_dataset{4},'a_dis');
a=read_from_file(LS_model_dataset{4},'a_dis');
yya = a ./ a_ref;
semilogx(f,yya(:,mode),'DisplayName', 'LS model');hold all

a_ref=read_from_file(const_model_dataset{4},'a_dis');
a=read_from_file(cigre_model_dataset{4},'a_dis');
yya = a ./ a_ref;
semilogx(f,yya(:,mode),'DisplayName', 'Cigre model');hold all

title(['mode # ' num2str(mode) ' - Attenuation constant, \rho = 100 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(1,2,2)

vel_ref=read_from_file(const_model_dataset{4},'vel_dis');
vel=read_from_file(LS_model_dataset{4},'vel_dis');
yyv = vel ./ vel_ref;
semilogx(f,yyv(:,mode),'DisplayName', 'LS model');hold all

vel_ref=read_from_file(const_model_dataset{4},'vel_dis');
vel=read_from_file(cigre_model_dataset{4},'vel_dis');
yyv = vel ./ vel_ref;
semilogx(f,yyv(:,mode),'DisplayName', 'Cigre model');hold all

title(['mode # ' num2str(mode) ' - Phase velocity, \rho = 100 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend

set_plot_params([1 2 2 1 1 4 1])
figure('Name',['FDCP___FDmodel___rho1000_mode' num2str(mode)])
subplot(1,2,1)

a_ref=read_from_file(const_model_dataset{5},'a_dis');
a=read_from_file(LS_model_dataset{5},'a_dis');
yya = a ./ a_ref;
semilogx(f,yya(:,mode),'DisplayName', 'LS model');hold all

a_ref=read_from_file(const_model_dataset{5},'a_dis');
a=read_from_file(cigre_model_dataset{5},'a_dis');
yya = a ./ a_ref;
semilogx(f,yya(:,mode),'DisplayName', 'Cigre model');hold all

title(['mode # ' num2str(mode) ' - Attenuation constant, \rho = 1000 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(1,2,2)

vel_ref=read_from_file(const_model_dataset{5},'vel_dis');
vel=read_from_file(LS_model_dataset{5},'vel_dis');
yyv = vel ./ vel_ref;
semilogx(f,yyv(:,mode),'DisplayName', 'LS model');hold all

vel_ref=read_from_file(const_model_dataset{5},'vel_dis');
vel=read_from_file(cigre_model_dataset{5},'vel_dis');
yyv = vel ./ vel_ref;
semilogx(f,yyv(:,mode),'DisplayName', 'Cigre model');hold all

title(['mode # ' num2str(mode) ' - Phase velocity, \rho = 1000 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend







set_plot_params([1 2 2 1 1 4 1])
figure('Name',['FDCP___FDmodel___rho10000_mode' num2str(mode)])
subplot(1,2,1)

a_ref=read_from_file(const_model_dataset{6},'a_dis');
a=read_from_file(LS_model_dataset{6},'a_dis');
yya = a ./ a_ref;
semilogx(f,yya(:,mode),'DisplayName', 'LS model');hold all

a_ref=read_from_file(const_model_dataset{6},'a_dis');
a=read_from_file(cigre_model_dataset{6},'a_dis');
yya = a ./ a_ref;
semilogx(f,yya(:,mode),'DisplayName', 'Cigre model');hold all

title(['mode # ' num2str(mode) ' - Attenuation constant, \rho = 5000 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\alpha_{norm} [p.u.]')
grid on
legend

subplot(1,2,2)

vel_ref=read_from_file(const_model_dataset{6},'vel_dis');
vel=read_from_file(LS_model_dataset{6},'vel_dis');
yyv = vel ./ vel_ref;
semilogx(f,yyv(:,mode),'DisplayName', 'LS model');hold all

vel_ref=read_from_file(const_model_dataset{6},'vel_dis');
vel=read_from_file(cigre_model_dataset{6},'vel_dis');
yyv = vel ./ vel_ref;
semilogx(f,yyv(:,mode),'DisplayName', 'Cigre model');hold all

title(['mode # ' num2str(mode) ' - Phase velocity, \rho = 5000 \Omega.m'])
axis tight
xlabel('Frequency [Hz]')
ylabel('\upsilon_{norm} [p.u.]')
grid on
legend