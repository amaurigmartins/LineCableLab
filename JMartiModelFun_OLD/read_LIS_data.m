% close all
% clear all
clc

ord=6;
for i=1:ord
    fname=['.\lisdata\mode' num2str(i) 'zc.txt'];
    data=load(fname);
    f=data(:,1);
    Zc=data(:,6);
    PHZC=data(:,7);
    [re,im] = pol2cart(PHZC*pi/180,Zc);
    mode(i).Zchvalue = complex(re,im);
    mode(i).Zchfreq=f;
    
    fname=['.\lisdata\mode' num2str(i) 'a1.txt'];
    data=load(fname);
    ff=data(:,1);
    vel = data(:,6)*1000;
    tau = data(:,7)*1e-3;
    A1=data(:,8);
    PHA1=data(:,9);
    [re,im] = pol2cart(PHA1*pi/180,A1);
    mode(i).A1value = complex(re,im);
    mode(i).A1freq = ff;
    mode(i).A1vel = vel;
    mode(i).A1tau = tau;
    
end

figure(9)
for i=1:ord
    semilogx(mode(i).Zchfreq,abs(mode(i).Zchvalue), 'DisplayName', ['mode #' num2str(i)]);hold all;
end
axis tight
xlabel('Frequency [Hz]')
ylabel('Z_{ch} [\Omega]')
grid on
legend;

figure(10)
for i=1:ord
    semilogx(mode(i).A1freq,abs(mode(i).A1value), 'DisplayName', ['mode #' num2str(i)]);hold all;
end
axis tight
xlabel('Frequency [Hz]')
ylabel('Propagation factor [unitless]')
grid on
legend;

% line_length = 1000;
% 
% % General flags
% spyplot_FLAG=0;
% msg_FLAG=true;
% ERR=1/100;
% stable_FLAG=0;
% 
% % Zc fit
% asymp_FLAG=2; %d!=0,e=0
% for m=1:ord
%     freq=mode(m).Zch.freq;
%     fun=mode(m).Zch.value;
%     [pol, res, infval, NORD, ffit, err]=vectfit_wrapper(fun,freq,ERR,asymp_FLAG,stable_FLAG,spyplot_FLAG,msg_FLAG);
%     fitZcATP(m).mode = m;
%     fitZcATP(m).NORD = NORD;
%     fitZcATP(m).zInf = infval;
%     fitZcATP(m).pol = pol;
%     fitZcATP(m).res = res;
%     figure(1);semilogx(freq,abs(ffit),'o','DisplayName', ['fit mode #' num2str(m)])
% end
% 
% % A1 fit
% spyplot_FLAG=0;
% asymp_FLAG=1; %d=e=0;
% for m=1:ord
%     tau=mode(m).A1.tau(end);
%     freq=mode(m).A1.freq;
%     fun = mode(m).A1.value.*exp(1i*2*pi.*freq.*tau);
%     [pol, res, infval, NORD, ffit, err]=vectfit_wrapper(fun,freq,ERR,asymp_FLAG,stable_FLAG,spyplot_FLAG,msg_FLAG);
%     fitA1ATP(m).mode = m;
%     fitA1ATP(m).NORD = NORD;
%     fitA1ATP(m).tauInf = tau;
%     fitA1ATP(m).pol = pol;
%     fitA1ATP(m).res = res;
%     figure(2);semilogx(freq,abs(ffit),'o','DisplayName', ['fit mode #' num2str(m)])
% end
% 
% fname=['.\lisdata\transf_matrix.mat'];
% load(fname);

% fname = 'E:\Users\Amauri\Documents\ATPdata\projects\Usp\pipefit.pch';
% fcontent = punchJMartiCard(ord, fitZcATP, fitA1ATP, Ti, fname);

