function [V]=freq_scan_sim(ord,freq,gamma_dis,Z_dis,Y_dis,Ti_dis,Ys_dis,Yr_dis,length_line,tag,NLTprnt)

if nargin == 10; NLTprnt=true;end

Vo1=ones(1,max(size(freq)));

gamma=gamma_dis;
Z=Z_dis;
Y=Y_dis;
Ti=Ti_dis;
Ys=Ys_dis;
Yr=Yr_dis;

Ybranch=zeros(max(size(freq)),(2*ord)^2); % (num_files x (2*ord)^2)   

for o=1:1:max(size(freq))
    Ybranch_dis=get_Ybranch_scaled(ord,length_line,Ti(o,:),Ys(o,:),Yr(o,:),gamma(o,:),Z(o,:)); % Function get_Ybranch - Calculates the admittance matrix at each frequency sample - (1 x (2*ord)^2)      
    Ybranch(o,:)=Ybranch_dis; % Stores each generated Ybranch_dis for each frequency sample - (num_files x (2*ord)^2)
end

%% Frequency Domain Calculation
V=calc_Vnode(ord,freq,Ybranch,Vo1); % Function calc_Vnode - Calculate node voltages fo S and R ends at FD  - (2*ord x max(size(f)))

%% Plot voltage spectrum
if NLTprnt
    Vmag=abs(V);

    for i=1:ord
        plottitle{i}=sprintf('Phase #%d - Sending terminal',i);
        plottitle{i+ord}=sprintf('Phase #%d - Receiving terminal',i);
    end

    for o=1:2*ord
        figure('Name', ['FreqScan' num2str(o) '_' tag])
        figure(o)
        semilogx(freq,Vmag(o,:),'LineWidth',2);
        xlabel('Frequency [Hz]')
        ylabel('Magnitude [pu]')
        title(plottitle{o});
        grid on
    end
end
