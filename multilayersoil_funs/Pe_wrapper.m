%% Wrapper for admittance function

function [Pe]=Pe_wrapper(EarthModel,h,d,r,e_g,m_g,sigma_g,t,f,con,flag,kx)

if nargin==11
    kx=0; %sets default to zero
end

% Constants
sig0=0;
eps0=8.8541878128e-12;
mu0=4*pi*1e-7;
w=2*pi*f;

e_g = [eps0; e_g];
m_g = [mu0; m_g];
sigma_g = [sig0; sigma_g];

if kx==0
    k_x= w.*0;
else
    k_x= w.*sqrt(m_g(kx).*(e_g(kx)-1i.*(sigma_g(kx)./w)));
end

num_layers=length(sigma_g);

for i=1:num_layers
    gamma(i)=sqrt(1i.*w.*m_g(i).*(sigma_g(i)+1i.*w.*e_g(i)));
    a{i}=@(lambda) sqrt(lambda.^2+gamma(i).^2+k_x.^2);
end


TOL=1e-3;

fun=matlabFunction(EarthModel.Y);
if isa(EarthModel.f,'sym')
    f_o=matlabFunction(EarthModel.f);
end

if isa(EarthModel.g,'sym')
    g_o=matlabFunction(EarthModel.g);
end

if isa(EarthModel.p,'sym')
    p_o=matlabFunction(EarthModel.p);
end

if isa(EarthModel.q,'sym')
    q_o=matlabFunction(EarthModel.q);
end
o=EarthModel.obs;

% construct integrand expression
if isa(EarthModel.f,'double')
    fun=@(lambda,hi,hj,yy) sum([0 ...
        fun(a{o}(lambda),...
        g_o(a{1}(lambda),a{2}(lambda),a{3}(lambda),gamma(o),hi,m_g(1),m_g(2),m_g(3),w,t),...
        gamma(o),hi,hj,lambda,m_g(o),w,...
        q_o(a{1}(lambda),a{2}(lambda),a{3}(lambda),gamma(1),gamma(2),gamma(3),hi,m_g(1),m_g(2),m_g(3),w,t),...
        yy)...
        ],'omitnan');
elseif isa(EarthModel.g,'double')
    fun=@(lambda,hi,hj,yy) sum([0 ...
        fun(a{o}(lambda),...
        f_o(a{1}(lambda),a{2}(lambda),a{3}(lambda),gamma(o),hi,m_g(1),m_g(2),m_g(3),w,t),...
        gamma(o),hi,hj,lambda,m_g(o),w,...
        p_o(a{1}(lambda),a{2}(lambda),a{3}(lambda),gamma(1),gamma(2),gamma(3),hi,m_g(1),m_g(2),m_g(3),w,t),...
        yy)...
        ],'omitnan');
else
    fun=@(lambda,hi,hj,yy) sum([0 ...
        fun(a{o}(lambda),...
        f_o(a{1}(lambda),a{2}(lambda),a{3}(lambda),gamma(o),hi,m_g(1),m_g(2),m_g(3),w,t),...
        g_o(a{1}(lambda),a{2}(lambda),a{3}(lambda),gamma(o),hi,m_g(1),m_g(2),m_g(3),w,t),...
        gamma(o),hi,hj,lambda,m_g(o),w,...
        p_o(a{1}(lambda),a{2}(lambda),a{3}(lambda),gamma(1),gamma(2),gamma(3),hi,m_g(1),m_g(2),m_g(3),w,t),...
        q_o(a{1}(lambda),a{2}(lambda),a{3}(lambda),gamma(1),gamma(2),gamma(3),hi,m_g(1),m_g(2),m_g(3),w,t),...
        yy)...
        ],'omitnan');
end
Pe=zeros(con,con);

if strcmp(flag,'mutual')
    for i=1:con
        for j=1:con
            if i~=j
                h1=h(1,i);
                h2=h(1,j);
                if abs(h2-h1)<TOL; h2=h2+TOL;end
                y=d(i,j);
                Qm=integral(@(lambda) fun(lambda,h1,h2,y),0,Inf,'ArrayValued',true);
                Pe(i,j)=1i*w*Qm;
            end
        end
    end
elseif strcmp(flag,'self')
    for k=1:1:con
        h1=h(1,k);
        h2=h(1,k)+TOL;
        y=r(k);
        Qm=integral(@(lambda) fun(lambda,h1,h2,y),0,Inf,'ArrayValued',true);
        Pe(k,k)=1i*w*Qm;
    end
end

% this ends the function. Remove it if you are dumb.
end