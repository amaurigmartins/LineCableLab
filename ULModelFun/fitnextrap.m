function [out_f, out_fun] = fitnextrap(f,in_fun,ord)

lastdec=log10(f(end));
new_freq=logspace(lastdec,lastdec+4,500);
new_freq_siz=length(f)+length(new_freq(2:end));
ff=[f; new_freq(2:end)'];
out_f=ff;
% ww = 2.*pi.*ff; % Frequency (rad/seg)
% ss = 1j*ww; % Complex Frequency
for o=1:ord
    fun=in_fun(:,o);
    exp_model = @(a,b,x) (a.*exp(-(x.*b)));
    obj_fun = @(params) norm(exp_model(params(1),params(2),f)-fun);
    sol = fminsearch(obj_fun, [1,0]);
    a_sol = sol(1);
    b_sol = sol(2);
    fun_extrap = exp_model(a_sol, b_sol, new_freq);
    %     mode(o).P=[fun; fun_extrap(2:end)'];
    out_fun(:,o)=[fun; fun_extrap(2:end)'];
          figure(10);semilogx(f,abs(fun),'ko');hold all;semilogx([f; new_freq(2:end)'],abs(out_fun));
    %     M=repmat(mode(o).Dj(:,:,freq_siz),[1 1 new_freq_siz-freq_siz-1]);
    %     mode(o).Dj(:,:,freq_siz+1:new_freq_siz)=repmat(mode(o).Dj(:,:,freq_siz),[1 1 new_freq_siz-freq_siz]);
end
TOL=1e-6;
out_fun=reshape([out_fun],[new_freq_siz,ord]);
out_fun(out_fun<=TOL)=0; %to prevent underflow issues

end