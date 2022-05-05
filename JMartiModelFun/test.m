dy=abs(fun(end)-fun(end-1));
df=abs(f(end)-f(end-1));
ff=f;
ffun=fun;

c=1;
while true
    newval=ffun(end)-dy;
    newf=ff(end)+df;
    
    dy=abs(ffun(end)-ffun(end-1))*10;
    df=abs(ff(end)-ff(end-1))*10;

    c=c+1;

    if newval <= .1
        break
    else
        ffun(end+1)=newval;
        ff(end+1)=newf;
    end
        
end

figure;semilogx(ff,abs(ffun));hold all; semilogx(f,abs(fun))

tol=2;
[P,Z,k] = Bode_process(ffun,ff,length(ff),tol);

as = k.*poly(Z); bs = poly(P); % Polynomials
[r,p,ks] = residue(as,bs) % Poles, residues and constant term

TF=isempty(ks);if(TF==1);ks=0;end
ffit = zeros(1,Ns)';
for k = 1:length(p)
    ffit = ffit + (r(k)./(s.' - p(k))).';
end
ffit = ffit + ks;

figure;semilogx(f,abs(ffit), 'o', 'DisplayName', ['fit mode #' num2str(m)])