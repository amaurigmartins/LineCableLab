function [] = plot_td_sim_fun(v,ord,sim_flag,data_t_sim,time_cl_brkr,tag)

supr_zero=find(data_t_sim==time_cl_brkr);

for i=1:ord
    plottitle{i}=sprintf('Phase #%d - Sending terminal',i);
    plottitle{i+ord}=sprintf('Phase #%d - Receiving terminal',i);
end

for o=1:2*ord
    v(1:supr_zero,o)=0;
    figure('Name', ['TDSim' num2str(o) '_' tag])
    plot(data_t_sim,v(:,o),'LineWidth',2);
    title(plottitle{o});
    xlabel('Time [s]')
    ylabel('Magnitude [pu]')
    grid on
    if sim_flag==4; xlim([0 50e-6]); end
end

end