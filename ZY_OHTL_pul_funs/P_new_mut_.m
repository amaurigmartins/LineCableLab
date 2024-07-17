function [Pg_mutual]=P_new_mut_(h,d,eps1,mu1,sigma1,f,con,kx,flag)


% THIS IS A MODIFIED VERSION THAT INCLUDES THE FORMULATION BASED ON THE MAGNETIC HERTZIAN


if nargin==8
    flag=0; % 0 for standard electric Hertzians, 1 for electric + magnetic
elseif nargin==7
    kx=0; %sets default to zero
    flag=0;
end

% Function for the Mutual Potential Coefficients with overhead-underground
% arrangement using the new formula

% Inputs
% h       : depth of conductor [m]
% d       : distances between conductors [m]
% e_g     : permittivity of earth [F/m]
% m_g     : permeability of earth [H/m]
% sigma_g : conductivity of earth [S/m]
% omega   : angular frequency [rad/s]

% Output
% out: Mutual earth potential coefficients between overhead/underground conductors [Ohm/m]

% Constants
sig0=0;
eps0=8.8541878128e-12;
mu0=4*pi*1e-7;
w=2*pi*f;

if flag==0
    if strcmp(kx,'k1')
        k_x=@(omega) omega.*sqrt(mu1.*(eps1-1i.*(sigma1./omega)));
    elseif strcmp(kx,'k0')
        k_x=@(omega) omega.*sqrt(mu0.*eps0);
    else
        k_x=@(omega) omega.*0;
    end

    gamma_0=@(omega) sqrt(1i.*omega.*mu0.*(sig0+1i.*omega.*eps0));
    gamma_1=@(omega) sqrt(1i.*omega.*mu1.*(sigma1+1i.*omega.*eps1));
    a_0=@(lambda,omega) sqrt(lambda.^2+gamma_0(omega).^2+k_x(omega).^2);
    a_1=@(lambda,omega) sqrt(lambda.^2+gamma_1(omega).^2+k_x(omega).^2);
elseif flag==1
    s1=(sig0+1i.*w.*eps0);
    s2=(sigma1+1i.*w.*eps1);
    mu2=mu1;
    mu1=mu0;
    gamma2= sqrt(1i.*w.*mu2.*s2);
    gamma1= sqrt(1i.*w.*mu1.*s1);
    a2=@(lambda) sqrt(lambda.^2+gamma2.^2);
    a1=@(lambda) sqrt(lambda.^2+gamma1.^2);
end

Pg_mutual=zeros(con,con);

for x=1:con
    for y=x+1:con
        if x~=y
            h1=h(1,x);
            h2=h(1,y);
            if flag==0
                if (h1 < 0 && h2 > 0) %Y10
                    yy=@(a0,a1,gamma0,gamma1,hi,hj,lambda,mu0,mu1,omega,y)(mu0.*mu1.*omega.*exp(a1.*hi-a0.*hj).*cos(lambda.*y).*(sign(hi)-1.0).*(a0.*mu0+a1.*mu1).*-5.0e-1i)./(pi.*(a0.*gamma1.^2.*mu0+a1.*gamma0.^2.*mu1).*(a0.*mu1+a1.*mu0));
                elseif (h1 > 0 && h2 < 0) %Y01
                    yy=@(a0,a1,gamma0,gamma1,hi,hj,lambda,mu0,mu1,omega,y)(mu0.*mu1.*omega.*exp(-a0.*hi+a1.*hj).*cos(lambda.*y).*(sign(hi)+1.0).*(a0.*mu0+a1.*mu1).*5.0e-1i)./(pi.*(a0.*gamma1.^2.*mu0+a1.*gamma0.^2.*mu1).*(a0.*mu1+a1.*mu0));
                else
                    continue
                end
                
                yfun=@(lambda,omega) sum([0 yy(a_0(lambda,omega),a_1(lambda,omega),gamma_0(omega),gamma_1(omega),h1,h2,lambda,mu0,mu1,omega,d(x,y))],'omitnan');
                yfun=@(lambda) yfun(lambda,w);
                Qm=integral(yfun,0,Inf,'ArrayValued',true);
                Pg_mutual(x,y)=1i*w*Qm;
            elseif flag==1
                % yy=@(lambda,hi,mu1,mu2,s2,a1,a2,y,hj)(mu1.*exp(-hi.*a1+a2.*hj).*cos(lambda.*y))./(s2.*pi.*(mu1.*a2+mu2.*a1));
                % yy=@(lambda,hi,mu1,mu2,s2,u1,u2,y,hj)(mu1.*exp(lambda.*y.*-1i-hi.*u1).*exp(u2.*hj))./(s2.*pi.*(mu1.*u2+mu2.*u1).*2.0);
                yy=@(L,h,m1,m2,s2,u1,u2,y,z)(m1.*exp(-h.*u1+u2.*z).*cos(L.*y).*2.0)./(s2.*pi.*(m1.*u2+m2.*u1));
                if (h1 < 0 && h2 > 0) %exploit reciprocity theorem and hope for the best
                    yfun=@(lambda) sum([0 yy(lambda,h2,mu1,mu2,s2,a1(lambda),a2(lambda),d(x,y),h1)],'omitnan');
                elseif (h1 > 0 && h2 < 0) 
                    yfun=@(lambda) sum([0 yy(lambda,h1,mu1,mu2,s2,a1(lambda),a2(lambda),d(x,y),h2)],'omitnan');
                else
                    continue
                end
            Qm=integral(yfun,-Inf,Inf,'ArrayValued',true);
            Pg_mutual(x,y)=1i*w*Qm;
            Pg_mutual(y,x)=Pg_mutual(x,y);


            end

            
            %             dij = sqrt((h1-h2)^2+d(x,y)^2);
            %             Dij = sqrt((h1+h2)^2+d(x,y)^2);
            %             Pg_mutual(x,y) = 1/(2*pi*eps0)*(log(Dij/dij));
        end
    end
end

end

