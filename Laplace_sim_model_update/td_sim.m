function [v]=td_sim(ord,f,freq,gamma_dis,Z_dis,Y_dis,Ti_dis,Ys_dis,Yr_dis,length,time_cl_brkr,samples,e,data_t_sim,time_sim,amp,sim_flag)


%% Build Circuit Struct - Fitting (4B)
%[gamma,Z,Y,Ti,Ys,Yr]=do_spline(ord,freq,f,gamma_dis,Z_dis,Y_dis,Ti_dis,Ys_dis,Yr_dis);
[gamma,Z,Y,Ti,Ys,Yr]=do_pchip(ord,freq,f,gamma_dis,Z_dis,Y_dis,Ti_dis,Ys_dis,Yr_dis);

Ybranch=zeros(max(size(f)),(2*ord)^2); % (num_files x (2*ord)^2)   

for o=1:1:max(size(f))
    % Ybranch_dis=get_Ybranch(ord,length,Ti(o,:),Ys(o,:),Yr(o,:),gamma(o,:),Z(o,:)); % Function get_Ybranch - Calculates the admittance matrix at each frequency sample - (1 x (2*ord)^2)
    Ybranch_dis=get_Ybranch_scaled(ord,length,Ti(o,:),Ys(o,:),Yr(o,:),gamma(o,:),Z(o,:)); % Function get_Ybranch - Calculates the admittance matrix at each frequency sample - (1 x (2*ord)^2)

    Ybranch(o,:)=Ybranch_dis; % Store all Ybranch_dis at each given frequency - (num_files x (2*ord)^2)
end

%% Voltage Sources and Frequency Domain Calculation (5)

if sim_flag == 1
    % a) Sinusoidal Voltage, flag=1
    [vo1,Vo1,c]=src_sin(amp,data_t_sim); % Function src_sin
    V=calc_Vnode(ord,f,Ybranch,Vo1); % Function calc_Vnode - Calculate node voltages fo S and R ends at FD  - (2*ord x max(size(f)))
elseif sim_flag == 2
    % b) Energization, flag=2
    [vo1,Vo1,c]=src_enrg(amp,time_cl_brkr,data_t_sim,e); % Functioen src_enrg
    V=calc_Vnode(ord,f,Ybranch,Vo1); % Function calc_Vnode - Calculate node voltages fo S and R ends at FD  - (2*ord x max(size(f)))
elseif sim_flag == 3
    % c) Step Response, flag=3
    [vo1,Vo1,c]=src_step(amp,time_cl_brkr,samples,e,f,data_t_sim,time_sim); % Functioen src_enrg
    V=calc_Vnode(ord,f,Ybranch,Vo1); % Function calc_Vnode - Calculate node voltages fo S and R ends at FD  - (2*ord x max(size(f)))
elseif sim_flag == 4
    % d) Double Exponential, flag=4
    [vo1,Vo1,c]=src_dexp(amp,time_cl_brkr,time_sim,samples,e,f,data_t_sim); % Function src_dexp with FFT or NLT!!
    V=calc_Vnode(ord,f,Ybranch,Vo1); % Function calc_Vnode - Calculate node voltages fo S and R ends at FD  - (2*ord x max(size(f)))
elseif sim_flag == 5
    % e) Custom Source Voltage, flag=5
    [vo1,Vo1,c]=src_custom_volt(time_sim,samples,e,data_t_sim);
    V=calc_Vnode(ord,f,Ybranch,Vo1); % Function calc_Vnode - Calculate node voltages fo S and R ends at FD  - (2*ord x max(size(f)))
elseif sim_flag == 6
    % f) 3ph Source, flag=6
    [vo1,vo2,vo3,Vo1,Vo2,Vo3,c]=src_3ph(amp,time_cl_brkr,data_t_sim,e,samples,f,time_sim); % Function src_3ph with FFT or NLT!!
    V=calc_Vnode_3ph(ord,f,Ybranch,Vo1,Vo2,Vo3);
    %V=calc_Vnode(ord,f,Ybranch,Vo1); % Function calc_Vnode - Calculate node voltages fo S and R ends at FD  - (2*ord x max(size(f)))
elseif sim_flag > 70
    % g) Lighting Current Source, flag=7
    [io1,Io1,c]=src_custom_curr(time_sim,samples,e,data_t_sim,sim_flag);
    V=calc_Vnode_curr(ord,f,Ybranch,Io1); % Function calc_Vnode - Calculate node voltages fo S and R ends at FD  - (2*ord x max(size(f)))
end

% spare) All frequencies
%[Vo1,c]=src_all_freq(data_t_sim,samples,time_sim,f,e); % Function src_all_freq % Function src_all_freq with FFT or NLT!!


%% Inverse FD - Time Domain Calculation (8)
v=tm_dmn_clc(V,f,e,data_t_sim,c,ord); % Function tm_dmn_clc - Calculate node voltages for S and R ends at TD - (2*max(size(f)) x 2*ord)=(samples x 2*ord)

% supr_zero=find(data_t_sim==time_cl_brkr);

% if NLTprnt
%     for i=1:ord
%         plottitle{i}=sprintf('Phase #%d - Sending terminal',i);
%         plottitle{i+ord}=sprintf('Phase #%d - Receiving terminal',i);
%     end
% 
%     for o=1:2*ord
%         v(1:supr_zero,o)=0;
%         figure('Name', ['TDSim' num2str(o) '_' tag])
%         plot(data_t_sim,v(:,o),'LineWidth',2);
%         title(plottitle{o});
%         xlabel('Time [s]')
%         ylabel('Magnitude [pu]')
%         grid on
%         if sim_flag==4; xlim([0 50e-6]); end;
%     end
% end
