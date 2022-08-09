figure(1)
subplot(3,1,1),plot(t(1:15001),vCsa(1:15001),t_B(1:15001),vCsa_B(1:15001),data_t_sim(1:15001),v(1:15001,1))
subplot(3,1,2),plot(t(1:15001),vCsb(1:15001),t_B(1:15001),vCsb_B(1:15001),data_t_sim(1:15001),v(1:15001,3))
subplot(3,1,3),plot(t(1:15001),vCsc(1:15001),t_B(1:15001),vCsc_B(1:15001),data_t_sim(1:15001),v(1:15001,5))

figure(2)
subplot(3,1,1),plot(t(1:15001),vSsa(1:15001),t_B(1:15001),vSsa_B(1:15001),data_t_sim(1:15001),v(1:15001,2))
subplot(3,1,2),plot(t(1:15001),vSsb(1:15001),t_B(1:15001),vSsb_B(1:15001),data_t_sim(1:15001),v(1:15001,4))
subplot(3,1,3),plot(t(1:15001),vSsc(1:15001),t_B(1:15001),vSsc_B(1:15001),data_t_sim(1:15001),v(1:15001,6))

figure(3)
subplot(3,1,1),plot(t(1:15001),vCra(1:15001),t_B(1:15001),vCra_B(1:15001),data_t_sim(1:15001),v(1:15001,7))
subplot(3,1,2),plot(t(1:15001),vCrb(1:15001),t_B(1:15001),vCrb_B(1:15001),data_t_sim(1:15001),v(1:15001,9))
subplot(3,1,3),plot(t(1:15001),vCrc(1:15001),t_B(1:15001),vCrc_B(1:15001),data_t_sim(1:15001),v(1:15001,11))

figure(4)
subplot(3,1,1),plot(t(1:15001),vSra(1:15001),t_B(1:15001),vSra_B(1:15001),data_t_sim(1:15001),v(1:15001,8))
subplot(3,1,2),plot(t(1:15001),vSrb(1:15001),t_B(1:15001),vSrb_B(1:15001),data_t_sim(1:15001),v(1:15001,10))
subplot(3,1,3),plot(t(1:15001),vSrc(1:15001),t_B(1:15001),vSrc_B(1:15001),data_t_sim(1:15001),v(1:15001,12))