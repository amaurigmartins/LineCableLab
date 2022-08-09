close all;

% figure(1)
% subplot(3,2,1);semilogx(freq,real(Ybranch_dis_tot(:,1)),transpose(f),real(Ybranch(:,1)));
% subplot(3,2,3);semilogx(freq,real(Ybranch_dis_tot(:,3)),transpose(f),real(Ybranch(:,3)));
% subplot(3,2,5);semilogx(freq,real(Ybranch_dis_tot(:,12)),transpose(f),real(Ybranch(:,12)));
% 
% 
% subplot(3,2,2);semilogx(freq,imag(Ybranch_dis_tot(:,1)),transpose(f),imag(Ybranch(:,1)));
% subplot(3,2,4);semilogx(freq,imag(Ybranch_dis_tot(:,3)),transpose(f),imag(Ybranch(:,3)));
% subplot(3,2,6);semilogx(freq,imag(Ybranch_dis_tot(:,12)),transpose(f),imag(Ybranch(:,12)));


%     subplot(3,2,1);semilogx(freq,real(Ybranch_dis_tot(:,1)),transpose(f),real(Ybranch(:,1)));
%     subplot(3,2,3);semilogx(freq,real(Ybranch_dis_tot(:,12)),transpose(f),real(Ybranch(:,12)));   
%     subplot(3,2,5);semilogx(freq,real(Ybranch_dis_tot(:,23)),transpose(f),real(Ybranch(:,23)));
%     subplot(3,2,2);semilogx(freq,imag(Ybranch_dis_tot(:,1)),transpose(f),imag(Ybranch(:,1)));
%     subplot(3,2,4);semilogx(freq,imag(Ybranch_dis_tot(:,12)),transpose(f),imag(Ybranch(:,12)));
%     subplot(3,2,6);semilogx(freq,imag(Ybranch_dis_tot(:,23)),transpose(f),imag(Ybranch(:,23)));
   
%     subplot(2,2,1);semilogx(freq,real(Ybranch_dis_tot(:,1)),transpose(f),real(Ybranch(:,1)));
%     subplot(2,2,3);semilogx(freq,real(Ybranch_dis_tot(:,8)),transpose(f),real(Ybranch(:,8)));
%     subplot(2,2,2);semilogx(freq,imag(Ybranch_dis_tot(:,1)),transpose(f),imag(Ybranch(:,1)));
%     subplot(2,2,4);semilogx(freq,imag(Ybranch_dis_tot(:,8)),transpose(f),imag(Ybranch(:,8)));

% i=1;
% for o=1:2:677
%    freq_pl(i)=freq(o);
%    Ybranch_dis_tot_pl(i,:)=Ybranch_dis_tot(o,:);
%    i=i+1;
% end
% 
% 
%     subplot(2,2,1);semilogx(transpose(f),real(Ybranch(:,1)));
%     subplot(2,2,3);semilogx(transpose(f),real(Ybranch(:,8)));
%     subplot(2,2,2);semilogx(transpose(f),imag(Ybranch(:,1)));
%     subplot(2,2,4);semilogx(transpose(f),imag(Ybranch(:,8)));

for i=1:1:(2*ord)^2
    semilogx(f,unwrap(radtodeg(angle(Ybranch(:,i)))))
    hold on
end
hold off