function [out_f, out_fun] = fitnextrap(f,in_fun,ord, method)


fx=f(end);
FUNLIM=0.1;
TOL=1e-6;
NDEC=4;

lastdec=log10(f(end));

if method=='expdecay'
    new_freq=logspace(lastdec,lastdec+NDEC,(lastdec+NDEC)*2);
    new_freq_siz=length(f)+length(new_freq(2:end));
    ff=[f; new_freq(2:end)'];
    for o=1:ord
        fun=abs(in_fun(:,o));
        exp_model = @(a,b,x) (a.*exp(-(x.*b)));
        obj_fun = @(params) norm(exp_model(params(1),params(2),f)-fun);
        sol = fminsearch(obj_fun, [1,0]);
        a_sol = sol(1);
        b_sol = sol(2);
        fun_extrap = exp_model(a_sol, b_sol, new_freq);
        out_fun(:,o)=[fun; fun_extrap(2:end)'];
                figure(10);semilogx(f,abs(fun),'ko');hold all;semilogx(ff,abs(out_fun));
    end
    out_fun=reshape([out_fun],[new_freq_siz,ord]);
    out_fun(out_fun<=TOL)=0; %to prevent underflow issues
elseif method=='decdecay'
    for o=1:ord
        
        fun=abs(in_fun(:,o));
        A=fun(end);
        for k=1:NDEC
            new_freq(k)=10^(lastdec+k);
            fun_extrap(k)=A*(fx/new_freq(k))^(NDEC*.1);
        end
        ff=[f; new_freq(1:end)'];
        new_freq_siz=length(ff);
        out_fun(:,o)=[fun; fun_extrap(1:end)'];
                figure(11);semilogx(f,abs(fun),'ko');hold all;semilogx(ff,abs(out_fun));
    end
    out_fun=reshape([out_fun],[new_freq_siz,ord]);
    
end
out_f=ff;

end