function [] = plotZ_fun_ct(f,ord,Z_pul,jobid)

allVars=fieldnames(Z_pul);
legendLabels={...
    'Carson';
    'Noda';
    'Deri';
    'Sunde';
    'Pettersson';
    'Wise';
    'Semlyen';
    'Alvarado-Betancourt';
    ['Papadopoulos' char(10) '(underground)'];
    ['Martins-Papadopoulos-Chrysochos ' char(10) '(mixed overhead-underground)'];
    ['Pollaczek (mixed' char(10) 'overhead-underground)'];
    'Kikuchi';
    ['Sunde 2-layered'];
    ['Xue (underground)'];
    };

% Plot results
for j=1:ord
    if j==1; plotLabel='Self';else;plotLabel='Mutual';end
    figure('Name', sprintf('%sImpZ1%d_%s',plotLabel,j,jobid))
    for i=1:length(allVars)
        subplot(2,1,1)
        if ~all(Z_pul.(allVars{i})==0)
            semilogx(f,squeeze(abs(Z_pul.(allVars{i})(1,j,:))),'LineWidth',2,'DisplayName',legendLabels{i});hold all
        end
    end
    legend('-DynamicLegend','Location','best');
    ylabel('Magnitude [\Omega/m]')
    xlabel('Frequency [Hz]')
    grid on
    title(sprintf('%s impedance - Z1%d',plotLabel,j))

    for i=1:length(allVars)
        subplot(2,1,2)
        if ~all(Z_pul.(allVars{i})==0)
            semilogx(f,rad2deg(unwrap(squeeze(angle(Z_pul.(allVars{i})(1,j,:))))),'LineWidth',2,'DisplayName',legendLabels{i});hold all
        end
    end
    ylabel('Angle [º]')
    xlabel('Frequency [Hz]')
    grid on
end

% 
% % Self impedance (Z11)
% % set_plot_params()
% figure('Name', ['SelfImpZ11_' jobid])
% o=1;
% subplot(2,1,1)
% if ~all(Ztot_Carson==0)
%     loglog(f,squeeze(abs(Ztot_Carson(1,1,:))),'LineWidth',2);hold all
%     lgd{o} = 'Carson';
%     o=o+1;
% end
% 
% if ~all(Ztot_Noda==0)
%     loglog(f,squeeze(abs(Ztot_Noda(1,1,:))),'LineWidth',2);hold all
%     lgd{o} = 'Noda';
%     o=o+1;
% end
% 
% if ~all(Ztot_Deri==0)
%     loglog(f,squeeze(abs(Ztot_Deri(1,1,:))),'LineWidth',2);hold all
%     lgd{o} = 'Deri';
%     o=o+1;
% end
% 
% if ~all(Ztot_AlDe==0)
%     loglog(f,squeeze(abs(Ztot_AlDe(1,1,:))),'LineWidth',2);hold all
%     lgd{o} = 'Alvarado-Betancourt';
%     o=o+1;
% end
% 
% if ~all(Ztot_Sunde==0)
%     loglog(f,squeeze(abs(Ztot_Sunde(1,1,:))),'LineWidth',2);hold all
%     lgd{o} = 'Sunde';
%     o=o+1;
% end
% 
% if ~all(Ztot_Pettersson==0)
%     loglog(f,squeeze(abs(Ztot_Pettersson(1,1,:))),'LineWidth',2);hold all
%     lgd{o} = 'Pettersson';
%     o=o+1;
% end
% 
% if ~all(Ztot_Semlyen==0)
%     loglog(f,squeeze(abs(Ztot_Semlyen(1,1,:))),'LineWidth',2);hold all
%     lgd{o} = 'Semlyen';
%     o=o+1;
% end
% 
% if ~all(Ztot_Wise==0)
%     loglog(f,squeeze(abs(Ztot_Wise(1,1,:))),'LineWidth',2);hold all
%     lgd{o} = 'Wise';
%     o=o+1;
% end
% 
% if ~all(Ztot_Kik==0)
%     loglog(f,squeeze(abs(Ztot_Kik(1,1,:))),'LineWidth',2);hold all
%     lgd{o} = 'Kikuchi';
%     o=o+1;
% end
% 
% if ~all(Ztot_under==0)
%     loglog(f,squeeze(abs(Ztot_under(1,1,:))),'LineWidth',2);hold all
%     lgd{o} = ['Papadopoulos' char(10) '(underground)'];
%     o=o+1;
% end
% 
% if ~all(Ztot_over_under==0)
%     loglog(f,squeeze(abs(Ztot_over_under(1,1,:))),'LineWidth',2);hold all
%     lgd{o} = ['New formulas (mixed' char(10) 'overhead-underground)'];
%     o=o+1;
% end
% 
% if ~all(Ztot_Sunde2La==0)
%     loglog(f,squeeze(abs(Ztot_Sunde2La(1,1,:))),'LineWidth',2);hold all
%     lgd{o} = ['Sunde 2-layered'];
%     o=o+1;
% end
% 
% if ~all(Ztot_Xue==0)
%     loglog(f,squeeze(abs(Ztot_Xue(1,1,:))),'LineWidth',2);hold all
%     lgd{o} = ['Xue (underground)'];
%     o=o+1;
% end
% 
% 
% % if ~all(Ztot_OvUndPol==0)
% %     loglog(f,squeeze(abs(Ztot_OvUndPol(1,1,:))),'LineWidth',2);hold all
% %     lgd{o} = ['Pollaczek (mixed' char(10) 'overhead-underground)'];
% %     o=o+1;
% % end
% 
% %loglog(f,squeeze(abs(Ztot_Carson(1,1,:))), ...
% %    f,squeeze(abs(Ztot_Noda(1,1,:))), ...
% %    f,squeeze(abs(Ztot_Deri(1,1,:))), ...
% %    f,squeeze(abs(Ztot_AlDe(1,1,:))), ...
% %    f,squeeze(abs(Ztot_Sunde(1,1,:))), ...
% %    f,squeeze(abs(Ztot_Pettersson(1,1,:))), ...
% %    f,squeeze(abs(Ztot_Semlyen(1,1,:))), ...
% %    f,squeeze(abs(Ztot_Wise(1,1,:))), ...
% %    f,squeeze(abs(Ztot_under(1,1,:))), ...
% %    f,squeeze(abs(Ztot_over_under(1,1,:))),'LineWidth',2)
% xlabel('Frequency [Hz]')
% ylabel('Magnitude [\Omega/m]')
% % legend('Carson','Noda','Deri','Alvarado-Betancourt','Sunde','Pettersson','Semlyen','Wise','Underground','Overhead-underground')
% legend(lgd);
% grid on
% title('Self impedance - Z11')
% 
% o=1;
% subplot(2,1,2)
% if ~all(Ztot_Carson==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Carson(1,1,:))))),'LineWidth',2);hold all
%     lgd{o} = 'Carson';
%     o=o+1;
% end
% 
% if ~all(Ztot_Noda==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Noda(1,1,:))))),'LineWidth',2);hold all
%     lgd{o} = 'Noda';
%     o=o+1;
% end
% 
% if ~all(Ztot_Deri==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Deri(1,1,:))))),'LineWidth',2);hold all
%     lgd{o} = 'Deri';
%     o=o+1;
% end
% 
% if ~all(Ztot_AlDe==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_AlDe(1,1,:))))),'LineWidth',2);hold all
%     lgd{o} = 'Alvarado-Betancourt';
%     o=o+1;
% end
% 
% if ~all(Ztot_Sunde==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Sunde(1,1,:))))),'LineWidth',2);hold all
%     lgd{o} = 'Sunde';
%     o=o+1;
% end
% 
% if ~all(Ztot_Pettersson==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Pettersson(1,1,:))))),'LineWidth',2);hold all
%     lgd{o} = 'Pettersson';
%     o=o+1;
% end
% 
% if ~all(Ztot_Semlyen==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Semlyen(1,1,:))))),'LineWidth',2);hold all
%     lgd{o} = 'Semlyen';
%     o=o+1;
% end
% 
% if ~all(Ztot_Wise==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Wise(1,1,:))))),'LineWidth',2);hold all
%     lgd{o} = 'Wise';
%     o=o+1;
% end
% 
% if ~all(Ztot_Kik==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Kik(1,1,:))))),'LineWidth',2);hold all
%     lgd{o} = 'Kikuchi';
%     o=o+1;
% end
% 
% if ~all(Ztot_under==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_under(1,1,:))))),'LineWidth',2);hold all
%     lgd{o} = ['Papadopoulos' char(10) '(underground)'];
%     o=o+1;
% end
% 
% if ~all(Ztot_over_under==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_over_under(1,1,:))))),'LineWidth',2);hold all
%     lgd{o} = ['New formulas (mixed' char(10) 'overhead-underground)'];
%     o=o+1;
% end
% 
% if ~all(Ztot_Sunde2La==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Sunde2La(1,1,:))))),'LineWidth',2);hold all
%     lgd{o} = ['Sunde 2-layered'];
%     o=o+1;
% end
% 
% if ~all(Ztot_Xue==0)
%     loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Xue(1,1,:))))),'LineWidth',2);hold all
%     lgd{o} = ['Xue (underground)'];
%     o=o+1;
% end
% 
% 
% % if ~all(Ztot_OvUndPol==0)
% %     loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_OvUndPol(1,1,:))))),'LineWidth',2);hold all
% %     lgd{o} = ['Kikuchi (mixed' char(10) 'overhead-underground)'];
% %     o=o+1;
% % end
% 
% %loglog(f,rad2deg(squeeze(angle(Ztot_Carson(1,1,:)))), ...
% %    f,rad2deg(squeeze(angle(Ztot_Noda(1,1,:)))), ...
% %    f,rad2deg(squeeze(angle(Ztot_Deri(1,1,:)))), ...
% %    f,rad2deg(squeeze(angle(Ztot_AlDe(1,1,:)))), ...
% %    f,rad2deg(squeeze(angle(Ztot_Sunde(1,1,:)))), ...
% %    f,rad2deg(squeeze(angle(Ztot_Pettersson(1,1,:)))), ...
% %    f,rad2deg(squeeze(angle(Ztot_Semlyen(1,1,:)))), ...
% %    f,rad2deg(squeeze(angle(Ztot_Wise(1,1,:)))), ...
% %    f,rad2deg(squeeze(angle(Ztot_under(1,1,:)))), ...
% %    f,rad2deg(squeeze(angle(Ztot_over_under(1,1,:)))),'LineWidth',2)
% xlabel('Frequency [Hz]')
% ylabel('Angle [deg]')
% % legend('Carson','Noda','Deri','Alvarado-Betancourt','Sunde','Pettersson','Semlyen','Wise','Underground','Overhead-underground')
% legend(lgd);
% grid on
% 
% % Mutual impedance (Z12)
% 
% number=ord;
% 
% counter=2;
% % set_plot_params()
% figure('Name', ['MutualImpZ12_' jobid])
% 
% 
% % Zm_pet=squeeze(abs(Ztot_Pettersson(1,number,:)));
% 
% while (number-1>0)
% 
%     o=1;
%     subplot(2,1,1)
%     if ~all(Ztot_Carson==0)
%         loglog(f,squeeze(abs(Ztot_Carson(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = 'Carson';
%         o=o+1;
%     end
% 
%     if ~all(Ztot_Noda==0)
%         loglog(f,squeeze(abs(Ztot_Noda(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = 'Noda';
%         o=o+1;
%     end
% 
%     if ~all(Ztot_Deri==0)
%         loglog(f,squeeze(abs(Ztot_Deri(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = 'Deri';
%         o=o+1;
%     end
% 
%     if ~all(Ztot_AlDe==0)
%         loglog(f,squeeze(abs(Ztot_AlDe(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = 'Alvarado-Betancourt';
%         o=o+1;
%     end
% 
%     if ~all(Ztot_Sunde==0)
%         loglog(f,squeeze(abs(Ztot_Sunde(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = 'Sunde';
%         o=o+1;
%     end
% 
%     if ~all(Ztot_Pettersson==0)
%         loglog(f,squeeze(abs(Ztot_Pettersson(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = 'Pettersson';
%         o=o+1;
%     end
% 
%     if ~all(Ztot_Semlyen==0)
%         loglog(f,squeeze(abs(Ztot_Semlyen(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = 'Semlyen';
%         o=o+1;
%     end
% 
%     if ~all(Ztot_Wise==0)
%         loglog(f,squeeze(abs(Ztot_Wise(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = 'Wise';
%         o=o+1;
%     end
% 
%     if ~all(Ztot_Kik==0)
%         loglog(f,squeeze(abs(Ztot_Kik(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = 'Kikuchi';
%         o=o+1;
%     end
% 
%     if ~all(Ztot_under==0)
%         loglog(f,squeeze(abs(Ztot_under(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = ['Papadopoulos' char(10) '(underground)'];
%         o=o+1;
%     end
% 
%     if ~all(Ztot_over_under==0)
%         loglog(f,squeeze(abs(Ztot_over_under(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = ['New formulas (mixed' char(10) 'overhead-underground)'];
%         o=o+1;
%     end
% 
%     if ~all(Ztot_Sunde2La==0)
%         loglog(f,squeeze(abs(Ztot_Sunde2La(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = ['Sunde 2-layered'];
%         o=o+1;
%     end
% 
%     if ~all(Ztot_OvUndPol==0)
%         loglog(f,squeeze(abs(Ztot_OvUndPol(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = ['Pollaczek (mixed' char(10) 'overhead-underground)'];
%         o=o+1;
%     end
% 
%     if ~all(Ztot_Xue==0)
%         loglog(f,squeeze(abs(Ztot_Xue(1,number,:))),'LineWidth',2);hold all
%         lgd{o} = ['Xue (underground)'];
%         o=o+1;
%     end
% 
% 
% 
% 
%     xlabel('Frequency [Hz]')
%     ylabel('Magnitude [\Omega/m]')
%     legend(lgd);
%     grid on
%     title(['Mutual impedance - Z1',num2str(number)])
% 
%     o=1;
%     subplot(2,1,2)
% 
%     if ~all(Ztot_Carson==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Carson(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = 'Carson';
%         o=o+1;
%     end
% 
%     if ~all(Ztot_Noda==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Noda(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = 'Noda';
%         o=o+1;
%     end
% 
%     if ~all(Ztot_Deri==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Deri(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = 'Deri';
%         o=o+1;
%     end
% 
%     if ~all(Ztot_AlDe==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_AlDe(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = 'Alvarado-Betancourt';
%         o=o+1;
%     end
% 
%     if ~all(Ztot_Sunde==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Sunde(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = 'Sunde';
%         o=o+1;
%     end
% 
%     if ~all(Ztot_Pettersson==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Pettersson(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = 'Pettersson';
%         o=o+1;
%     end
% 
%     if ~all(Ztot_Semlyen==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Semlyen(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = 'Semlyen';
%         o=o+1;
%     end
% 
%     if ~all(Ztot_Wise==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Wise(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = 'Wise';
%         o=o+1;
%     end
% 
%     if ~all(Ztot_Kik==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Kik(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = 'Kikuchi';
%         o=o+1;
%     end
% 
%     if ~all(Ztot_under==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_under(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = ['Papadopoulos' char(10) '(underground)'];
%         o=o+1;
%     end
% 
%     if ~all(Ztot_over_under==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_over_under(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = ['New formulas (mixed' char(10) 'overhead-underground)'];
%         o=o+1;
%     end
% 
%     if ~all(Ztot_Sunde2La==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Sunde2La(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = ['Sunde 2-layered'];
%         o=o+1;
%     end
% 
% 
%     if ~all(Ztot_OvUndPol==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_OvUndPol(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = ['Pollaczek (mixed' char(10) 'overhead-underground)'];
%         o=o+1;
%     end
% 
%     if ~all(Ztot_Xue==0)
%         loglog(f,rad2deg(unwrap(squeeze(angle(Ztot_Xue(1,number,:))))),'LineWidth',2);hold all
%         lgd{o} = ['Xue (underground)'];
%         o=o+1;
%     end
% 
%     xlabel('Frequency [Hz]')
%     ylabel('Angle [deg]')
%     legend(lgd);
%     grid on
%     number=number-1;
%     counter=counter+1;
%     if counter<=ord
%         figure('Name', ['MutualImpZ1' num2str(number) '_' jobid])
%     end
% 
% end
end