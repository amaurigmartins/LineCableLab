function [] = plotY_fun_ct(f,ord,Y_pul,jobid)

allVars=fieldnames(Y_pul);
legendLabels={...
    'Images';
    'Pettersson';
    'Wise';
    ['Papadopoulos' char(10) '(underground)'];
    ['Martins-Papadopoulos-Chrysochos ' char(10) '(mixed overhead-underground)'];
    'Kikuchi';
    ['Sunde 2-layered'];
    ['Xue (underground)'];
    ['Classical TL ' char(10) '(overhead-underground neglected)'];
    };

% Plot results
for j=1:ord
    if j==1; plotLabel='Self';else;plotLabel='Mutual';end
    figure('Name', sprintf('%sAdmY1%d_%s',plotLabel,j,jobid))
    for i=1:length(allVars)
        subplot(2,1,1)
        if ~all(Y_pul.(allVars{i})==0)
            semilogx(f,squeeze(abs(Y_pul.(allVars{i})(1,j,:))),'LineWidth',2,'DisplayName',legendLabels{i});hold all
        end
    end
    legend('-DynamicLegend','Location','best');
    ylabel('Magnitude [Sm]')
    xlabel('Frequency [Hz]')
    grid on
    title(sprintf('%s admittance - Y1%d',plotLabel,j))

    for i=1:length(allVars)
        subplot(2,1,2)
        if ~all(Y_pul.(allVars{i})==0)
            semilogx(f,rad2deg(unwrap(squeeze(angle(Y_pul.(allVars{i})(1,j,:))))),'LineWidth',2,'DisplayName',legendLabels{i});hold all
        end
    end
    ylabel('Angle [º]')
    xlabel('Frequency [Hz]')
    grid on
end

% % Plot results
% % Self admittance (Y11)
% % set_plot_params()
% figure('Name', ['SelfAdmY11_' jobid])
% o=1;
% 
% subplot(2,1,1)
% con=2;
% if ~all(Ytot_Imag==0)
%     loglog(f,squeeze(abs(Ytot_Imag(con,con,:))),'LineWidth',2);hold all
%     lgd{o} = 'Image';
%     o=o+1;
% end
% 
% if ~all(Ytot_Pet==0)
%     loglog(f,squeeze(abs(Ytot_Pet(con,con,:))),'LineWidth',2);hold all
%     lgd{o} = 'Pettersson';
%     o=o+1;
% end
% 
% if ~all(Ytot_Wise==0)
%     loglog(f,squeeze(abs(Ytot_Wise(con,con,:))),'LineWidth',2);hold all
%     lgd{o} = 'Wise';
%     o=o+1;
% end
% 
% if ~all(Ytot_Kik==0)
%     loglog(f,squeeze(abs(Ytot_Kik(con,con,:))),'LineWidth',2);hold all
%     lgd{o} = 'Kikuchi';
%     o=o+1;
% end
% 
% if ~all(Ytot_under==0)
%     loglog(f,squeeze(abs(Ytot_under(con,con,:))),'LineWidth',2);hold all
%     lgd{o} = ['Papadopoulos' char(10) '(underground)'];
%     o=o+1;
% end
% 
% if ~all(Ytot_over_under==0)
%     loglog(f,squeeze(abs(Ytot_over_under(con,con,:))),'LineWidth',2);hold all
%     lgd{o} = ['New formulas (mixed' char(10) 'overhead-underground)'];
%     o=o+1;
% end
% 
% if ~all(Ytot_Sunde2La==0)
%     loglog(f,squeeze(abs(Ytot_Sunde2La(con,con,:))),'LineWidth',2);hold all
%     lgd{o} = ['Sunde 2-layered'];
%     o=o+1;
% end
% 
% if ~all(Ytot_Xue==0)
%     loglog(f,squeeze(abs(Ytot_Xue(con,con,:))),'LineWidth',2);hold all
%     lgd{o} = ['Xue (underground)'];
%     o=o+1;
% end
% 
% xlabel('Frequency [Hz]')
% ylabel('Magnitude [Sm]')
% legend(lgd)
% grid on
% title(sprintf('Self admittance - Y%d%d',con,con))
% 
% subplot(2,1,2)
% o=1;
% 
% if ~all(Ytot_Imag==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Imag(con,con,:))))),'LineWidth',2);hold all
%     lgd{o} = 'Image';
%     o=o+1;
% end
% 
% if ~all(Ytot_Pet==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Pet(con,con,:))))),'LineWidth',2);hold all
%     lgd{o} = 'Pettersson';
%     o=o+1;
% end
% 
% if ~all(Ytot_Wise==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Wise(con,con,:))))),'LineWidth',2);hold all
%     lgd{o} = 'Wise';
%     o=o+1;
% end
% 
% if ~all(Ytot_Kik==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Kik(con,con,:))))),'LineWidth',2);hold all
%     lgd{o} = 'Kikuchi';
%     o=o+1;
% end
% 
% if ~all(Ytot_under==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_under(con,con,:))))),'LineWidth',2);hold all
%     lgd{o} = ['Papadopoulos' char(10) '(underground)'];
%     o=o+1;
% end
% 
% if ~all(Ytot_over_under==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_over_under(con,con,:))))),'LineWidth',2);hold all
%     lgd{o} = ['New formulas (mixed' char(10) 'overhead-underground)'];
%     o=o+1;
% end
% 
% if ~all(Ytot_Sunde2La==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Sunde2La(con,con,:))))),'LineWidth',2);hold all
%     lgd{o} = ['Sunde 2-layered'];
%     o=o+1;
% end
% 
% if ~all(Ytot_Xue==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Xue(con,con,:))))),'LineWidth',2);hold all
%     lgd{o} = ['Xue (underground)'];
%     o=o+1;
% end
% 
% xlabel('Frequency [Hz]')
% ylabel('Angle [deg]')
% legend(lgd)
% grid on
% 
% % Mutual admittance (Y12)
% 
% number=ord;
% 
% % Ym_pet=squeeze(abs(Ytot_Pet(1,number,:)));
% 
% counter=ord+2;
% % set_plot_params()
% figure('Name', ['MutualAdmY12_' jobid])
% 
% while (number-1>0)
% 
%     subplot(2,1,1)
%     o=1;
% 
%     if ~all(Ytot_Imag==0)
%         loglog(f,squeeze(abs(Ytot_Imag(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = 'Image';
%         o=o+1;
%     end
% 
%     if ~all(Ytot_Pet==0)
%         loglog(f,squeeze(abs(Ytot_Pet(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = 'Pettersson';
%         o=o+1;
%     end
% 
%     if ~all(Ytot_Wise==0)
%         loglog(f,squeeze(abs(Ytot_Wise(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = 'Wise';
%         o=o+1;
%     end
% 
%     if ~all(Ytot_Kik==0)
%         loglog(f,squeeze(abs(Ytot_Kik(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = 'Kikuchi';
%         o=o+1;
%     end
% 
% 
%     if ~all(Ytot_under==0)
%         loglog(f,squeeze(abs(Ytot_under(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = ['Papadopoulos' char(10) '(underground)'];
%         o=o+1;
%     end
% 
%     if ~all(Ytot_over_under==0)
%         loglog(f,squeeze(abs(Ytot_over_under(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = ['New formulas (mixed' char(10) 'overhead-underground)'];
%         o=o+1;
%     end
% 
%     if ~all(Ytot_Sunde2La==0)
%         loglog(f,squeeze(abs(Ytot_Sunde2La(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = ['Sunde 2-layered'];
%         o=o+1;
%     end
% 
%     if ~all(Ytot_Xue==0)
%         loglog(f,squeeze(abs(Ytot_Xue(1,number,:))),'LineWidth',2);hold all
%         % lgd{o} = ['Xue (underground)'];
%         o=o+1;
%     end
% 
% 
%     xlabel('Frequency [Hz]')
%     ylabel('Magnitude [Sm]')
%     legend(lgd)
%     grid on
%     title(['Mutual admittance - Y1',num2str(number)])
% 
%     subplot(2,1,2)
%     o=1;
% 
%     if ~all(Ytot_Imag==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Imag(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = 'Image';
%         o=o+1;
%     end
% 
%     if ~all(Ytot_Pet==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Pet(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = 'Pettersson';
%         o=o+1;
%     end
% 
%     if ~all(Ytot_Wise==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Wise(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = 'Wise';
%         o=o+1;
%     end
% 
%     if ~all(Ytot_Kik==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Kik(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = 'Kikuchi';
%         o=o+1;
%     end
% 
% 
%     if ~all(Ytot_under==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_under(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = ['Papadopoulos' char(10) '(underground)'];
%         o=o+1;
%     end
% 
%     if ~all(Ytot_over_under==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_over_under(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = ['New formulas (mixed' char(10) 'overhead-underground)'];
%         o=o+1;
%     end
% 
%     if ~all(Ytot_Sunde2La==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Sunde2La(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = ['Sunde 2-layered'];
%         o=o+1;
%     end
% 
%     if ~all(Ytot_Xue==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ytot_Xue(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = ['Xue (underground)'];
%         o=o+1;
%     end
% 
%     xlabel('Frequency [Hz]')
%     ylabel('Angle [deg]')
%     legend(lgd)
%     grid on
%     number=number-1;
%     counter=counter+1;
%     if counter<=2*ord
%         figure('Name', ['MutualAdmY1' num2str(number) '_' jobid])
%     end

end
