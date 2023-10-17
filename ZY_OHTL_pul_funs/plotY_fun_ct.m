function [] = plotY_fun_ct(f,ord,Ytot_Imag,Ytot_Pet,Ytot_Wise, Ytot_under, Ytot_over_under,Ytot_Kik,Ytot_Sunde2La,jobid)

%% Plot results
% Self admittance (Y11)
set_plot_params()
figure('Name', ['SelfAdmY11_' jobid])
o=1;

subplot(2,1,1)

if ~all(Ytot_Imag==0)
    loglog(f,squeeze(abs(Ytot_Imag(1,1,:))),'LineWidth',2);hold all
    lgd{o} = 'Image';
    o=o+1;
end

if ~all(Ytot_Pet==0)
    loglog(f,squeeze(abs(Ytot_Pet(1,1,:))),'LineWidth',2);hold all
    lgd{o} = 'Pettersson';
    o=o+1;
end

if ~all(Ytot_Wise==0)
    loglog(f,squeeze(abs(Ytot_Wise(1,1,:))),'LineWidth',2);hold all
    lgd{o} = 'Wise';
    o=o+1;
end

if ~all(Ytot_Kik==0)
    loglog(f,squeeze(abs(Ytot_Kik(1,1,:))),'LineWidth',2);hold all
    lgd{o} = 'Kikuchi';
    o=o+1;
end

if ~all(Ytot_under==0)
    loglog(f,squeeze(abs(Ytot_under(1,1,:))),'LineWidth',2);hold all
    lgd{o} = ['Papadopoulos' char(10) '(underground)'];
    o=o+1;
end

if ~all(Ytot_over_under==0)
    loglog(f,squeeze(abs(Ytot_over_under(1,1,:))),'LineWidth',2);hold all
    lgd{o} = ['New formulas (mixed' char(10) 'overhead-underground)'];
    o=o+1;
end

if ~all(Ytot_Sunde2La==0)
    loglog(f,squeeze(abs(Ytot_Sunde2La(1,1,:))),'LineWidth',2);hold all
    lgd{o} = ['Sunde 2-layered'];
    o=o+1;
end

xlabel('Frequency [Hz]')
ylabel('Magnitude [Sm]')
legend(lgd)
grid on
title('Self admittance - Y11')

subplot(2,1,2)
o=1;

if ~all(Ytot_Imag==0)
    loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Imag(1,1,:))))),'LineWidth',2);hold all
    lgd{o} = 'Image';
    o=o+1;
end

if ~all(Ytot_Pet==0)
    loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Pet(1,1,:))))),'LineWidth',2);hold all
    lgd{o} = 'Pettersson';
    o=o+1;
end

if ~all(Ytot_Wise==0)
    loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Wise(1,1,:))))),'LineWidth',2);hold all
    lgd{o} = 'Wise';
    o=o+1;
end

if ~all(Ytot_Kik==0)
    loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Kik(1,1,:))))),'LineWidth',2);hold all
    lgd{o} = 'Kikuchi';
    o=o+1;
end

if ~all(Ytot_under==0)
    loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_under(1,1,:))))),'LineWidth',2);hold all
    lgd{o} = ['Papadopoulos' char(10) '(underground)'];
    o=o+1;
end

if ~all(Ytot_over_under==0)
    loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_over_under(1,1,:))))),'LineWidth',2);hold all
    lgd{o} = ['New formulas (mixed' char(10) 'overhead-underground)'];
    o=o+1;
end

if ~all(Ytot_Sunde2La==0)
    loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Sunde2La(1,1,:))))),'LineWidth',2);hold all
    lgd{o} = ['Sunde 2-layered'];
    o=o+1;
end

xlabel('Frequency [Hz]')
ylabel('Angle [deg]')
legend(lgd)
grid on

% Mutual admittance (Y12)

number=ord;

% Ym_pet=squeeze(abs(Ytot_Pet(1,number,:)));

counter=ord+2;
set_plot_params()
figure('Name', ['MutualAdmY12_' jobid])

while (number-1>0)

    subplot(2,1,1)
    o=1;

    if ~all(Ytot_Imag==0)
        loglog(f,squeeze(abs(Ytot_Imag(1,number,:))),'LineWidth',2);hold all
        lgd{o} = 'Image';
        o=o+1;
    end

    if ~all(Ytot_Pet==0)
        loglog(f,squeeze(abs(Ytot_Pet(1,number,:))),'LineWidth',2);hold all
        lgd{o} = 'Pettersson';
        o=o+1;
    end

    if ~all(Ytot_Wise==0)
        loglog(f,squeeze(abs(Ytot_Wise(1,number,:))),'LineWidth',2);hold all
        lgd{o} = 'Wise';
        o=o+1;
    end

    if ~all(Ytot_Kik==0)
        loglog(f,squeeze(abs(Ytot_Kik(1,number,:))),'LineWidth',2);hold all
        lgd{o} = 'Kikuchi';
        o=o+1;
    end


    if ~all(Ytot_under==0)
        loglog(f,squeeze(abs(Ytot_under(1,number,:))),'LineWidth',2);hold all
        lgd{o} = ['Papadopoulos' char(10) '(underground)'];
        o=o+1;
    end

    if ~all(Ytot_over_under==0)
        loglog(f,squeeze(abs(Ytot_over_under(1,number,:))),'LineWidth',2);hold all
        lgd{o} = ['New formulas (mixed' char(10) 'overhead-underground)'];
        o=o+1;
    end

    if ~all(Ytot_Sunde2La==0)
        loglog(f,squeeze(abs(Ytot_Sunde2La(1,number,:))),'LineWidth',2);hold all
        lgd{o} = ['Sunde 2-layered'];
        o=o+1;
    end


    xlabel('Frequency [Hz]')
    ylabel('Magnitude [Sm]')
    legend(lgd)
    grid on
    title(['Mutual admittance - Y1',num2str(number)])

    subplot(2,1,2)
    o=1;

    if ~all(Ytot_Imag==0)
        loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Imag(1,number,:))))),'LineWidth',2);hold all
        lgd{o} = 'Image';
        o=o+1;
    end

    if ~all(Ytot_Pet==0)
        loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Pet(1,number,:))))),'LineWidth',2);hold all
        lgd{o} = 'Pettersson';
        o=o+1;
    end

    if ~all(Ytot_Wise==0)
        loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Wise(1,number,:))))),'LineWidth',2);hold all
        lgd{o} = 'Wise';
        o=o+1;
    end

    if ~all(Ytot_Kik==0)
        loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Kik(1,number,:))))),'LineWidth',2);hold all
        lgd{o} = 'Kikuchi';
        o=o+1;
    end


    if ~all(Ytot_under==0)
        loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_under(1,number,:))))),'LineWidth',2);hold all
        lgd{o} = ['Papadopoulos' char(10) '(underground)'];
        o=o+1;
    end

    if ~all(Ytot_over_under==0)
        loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_over_under(1,number,:))))),'LineWidth',2);hold all
        lgd{o} = ['New formulas (mixed' char(10) 'overhead-underground)'];
        o=o+1;
    end

    if ~all(Ytot_Sunde2La==0)
        loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Sunde2La(1,number,:))))),'LineWidth',2);hold all
        lgd{o} = ['Sunde 2-layered'];
        o=o+1;
    end

    xlabel('Frequency [Hz]')
    ylabel('Angle [deg]')
    legend(lgd)
    grid on
    number=number-1;
    counter=counter+1;
    if counter<=2*ord
        figure('Name', ['MutualAdmY1' num2str(number) '_' jobid])
    end

end
