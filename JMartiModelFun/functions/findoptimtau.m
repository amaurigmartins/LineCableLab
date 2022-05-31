function [tau_a, tau_b, tau_opt] = findoptimtau(f,g_dis,line_length,ord)
freq_siz=length(f);

vel(1:ord)=(2*pi*f(freq_siz))./imag(g_dis(freq_siz,:));

for o=1:ord
    H(:,o)=exp(-g_dis(:,o).*line_length);
end

[ff, H] = fitnextrap(f,H,ord);

j=freq_siz;
w=2*pi*ff;
Ns=length(w);

MINFUN=.1;

for o=1:ord
    absH=abs(H(:,o));

    %First term in (9): Accurate_transmission_line_modeling_through_optima.pdf
    phase1=(pi/2)*log((absH(j+1)/absH(j-1)))/(log(w(j+1)/w(j-1)));
    %Second term in (9):
    phase2=0;
    term2=log((absH(j+1)/absH(j-1))) /(log(w(j+1)/w(j-1)));
    for k=2:Ns-1
        term1=log(absH(k+1)/absH(k-1)) /(log(w(k+1)/w(k-1)));
        if k~=j
            phase2=phase2+(abs(term1)-abs(term2))*log(coth(abs(log(w(k)/w(j)))/2))*log(w(k+1)/w(k));
        end
    end

    phase2=phase2/pi;
    ph=(phase1+phase2);
    if ph<0
        ph=ph+(2*pi);
    end
    phase_min(:,o)=ph; %Phase angle [rad]
end

tau_a(1:ord)=line_length./vel;
tau_b(1:ord)=line_length./vel + phase_min/(w(j));
options = optimset('Display','none', 'TolX',1e-6);

for o=1:ord
    tau_opt(o)=fminbnd(@(x)fcalc(H(:,o),w,x,10),tau_a(o),tau_b(o),options);
end

end

function out = fcalc(H, w, tau, N)
s=1i.*w;
fun = H.*exp(s.*tau);
Ns=length(s);
poles=logspace(s(1),s(end),N);
weight=ones(1,Ns);
opts.relax=1;
opts.stable=1;
opts.asymp=1;
opts.spy1=0;
opts.spy2=0;
opts.cmplx_ss=1;

if size(s,1) > 1
    s=s';
end

if size(fun,1) > 1
    fun=fun';
end

for i=1:N
    [~,poles,rmserr,~,~]=vectfit3(fun,s,poles,weight,opts);
end

out=rmserr;

end
