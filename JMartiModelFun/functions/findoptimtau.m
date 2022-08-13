function [rmserr, pol, res, ks, NORD, ffit, tau_mps, tau_opt, tau] = findoptimtau(f,vel,H,line_length,ERR,opts)
j=length(f)-1;

% [ff, H] = fitnextrap(f,H);
ff=f;

w=2*pi*ff;
Ns=length(w);

absH=abs(H);

w=2*pi*ff;

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
phase_min=ph; %Phase angle [rad]
tau_mps=(line_length/vel(j)) + phase_min/(w(j));


tau_a=line_length/3e8;
tau_b=(line_length/vel(j));


tmin=.8*tau_a;
tmax=1.2*tau_b;
tau=tau_b;

% DEBUGME
% dt=abs(tmin-tmax)/100;
% ttau=tmin:dt:tmax
% % opts.NORD=20;
% % opts.Niter=20;
% for i=1:length(ttau)
%     [rmserr, pol, res, ks, NORD, ffit]=fcalc(H, f, ttau(i), ERR, opts);pol
%     err_tau(i)=rmserr;
% end
% figure;plot(ttau,err_tau)

options = optimset('Display','none', 'TolX',1e-6);
tau_opt=fminbnd(@(x)fcalc(H,ff,x,ERR,opts),tmin,tmax,options);
opts.NORD=20;
opts.Niter=20;


while 1

    tau_opt=fminbnd(@(x)fcalc(H,ff,x,ERR,opts),tmin,tmax,options);
%     tau_opt=tau_b;
    [rmserr, pol, res, ks, NORD, ffit]=fcalc(H, ff, tau_opt, ERR, opts);

    if max(imag(pol))  ~= 0
        opts.NORD = opts.NORD-1;
    else
        break
    end
end

end


function [out, pol, res, ks, NORD, ffit] = fcalc(H, f, tau, ERR, opts)
fun = H.*exp(1i.*2.*pi.*f.*tau);
[pol, res, ks, NORD, ffit, err] = vectfit_wrapper(fun,f,ERR,opts);
out=err;
% pol
end
