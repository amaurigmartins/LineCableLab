function [Z]=Z_femm_slf_mut(Geom,soilFD,k,freq,con,basename)

Z = zeros(con,con);

try; closefemm; end

for i=1:con
    %% Open a FEMM instance
    fclose('all');
    
    fname=sprintf('%s_c%d.fem',basename,i);
    global HandleToFEMM;
    openfemm(1);
    hand1=HandleToFEMM;
    % create new document
    DOCTYPE=0; %0 for magnetics, 1 for electrostatics, 3 for current flow
    newdocument(DOCTYPE);
    main_minimize;
    
    %% Set problem solver
    w=2*pi*freq;
    minangle=33.8;
    depth=1;
    acsolver=1; %0 for successive approximation, 1 for Newton.

    mi_smartmesh(1)
    
    mi_probdef(freq,'meters','planar',1e-12,depth,minangle,acsolver);
    %mi_setgrid(1,'cart')
    mu0=4*pi*1e-7;
    
    
    %% Define energizations
    Isrc=1;
    for j=1:con
        if j==i
            mi_addcircprop(sprintf('cond%d_source',j), Isrc, 0);
        else
            mi_addcircprop(sprintf('cond%d_target',j), 0, 0);
        end
    end
    
    %% Define materials
    % Air and earth layers
    mi_addmaterial('air',1,1,0,0,0,0,0,1,0,0,0,0,0) % air layer
    if isfield(soilFD,'layer')
        % will handle this when I handle this
    else
        m_g=soilFD.m_g/mu0;
        % eps_g=soilFD.erg_total(k);
        sig_g=soilFD.sigma_g_total(k)/1e6;
        mi_addmaterial('layer1', m_g, m_g, 0, 0, sig_g, 0, 0, 1, 0, 0, 0, 0, 0) % earth layer
    end
    
    % Conductors and insulation
    for j=1:con
        m_c=Geom(j,7);
        sig_c=1/Geom(j,6)/1e6;
        mi_addmaterial(sprintf('cond%d_core',j), m_c, m_c, 0, 0, sig_c, 0, 0, 1, 0, 0, 0, 0, 0) % conductor core
        if ~isnan(Geom(j,8))
            m_i=Geom(j,9);
            sig_i=0;
            %             eps_i=Geom(j,10);
            mi_addmaterial(sprintf('cond%d_insu',j), m_i, m_i, 0, 0, sig_i, 0, 0, 1, 0, 0, 0, 0, 0) % conductor insulation
        end
    end
    
    %% Position objects in the 2D grid space
    L=50e3; % domain size
    
    % Air layer
    mi_addblocklabel(0,L/2);
    mi_selectlabel(0,L/2);
    mi_setblockprop('air',1,1,0,0,0,1);
    mi_clearselected;
    
    % Earth layers
    if isfield(soilFD,'layer')
        % will handle this when I handle this
    else
        t=[0 -L]; % soil layer interfaces
        mi_addblocklabel(0,(t(end)-t(1))/2);
        mi_selectlabel(0,(t(end)-t(1))/2);
        mi_setblockprop('layer1',1,0,0,0,0,1);
        mi_clearselected;
    end
    
    % Conductors
    for j=1:con
        x_c = Geom(j,2);
        y_c = Geom(j,3);
        r_in = Geom(j,4);
        r_ext = Geom(j,5);
        r_ins = Geom(j,8);
        if r_in > 0
            %draw inner circle
            mi_drawarc(x_c-r_in,y_c,x_c+r_in,y_c,180,3)
            mi_drawarc(x_c+r_in,y_c,x_c-r_in,y_c,180,3)
            mi_addblocklabel(x_c,y_c);
            mi_selectlabel(x_c,y_c);
            mi_setblockprop('air',1,1,0,0,0,1);
        end
        
        % draw outer circle
        mi_drawarc(x_c-r_ext,y_c,x_c+r_ext,y_c,180,3)
        mi_drawarc(x_c+r_ext,y_c,x_c-r_ext,y_c,180,3)
        mi_addblocklabel(x_c+(r_in+r_ext)/2,y_c);
        mi_selectlabel(x_c+(r_in+r_ext)/2,y_c);
        if j==i
            mi_setblockprop(sprintf('cond%d_core',j),0,0.001,sprintf('cond%d_source',j),0,0,1)
        else
            mi_setblockprop(sprintf('cond%d_core',j),0,0.001,sprintf('cond%d_target',j),0,0,1)
        end
        mi_clearselected
        
        % draw insulation
        if ~isnan(r_ins)
            mi_drawarc(x_c-r_ins,y_c,x_c+r_ins,y_c,180,3)
            mi_drawarc(x_c+r_ins,y_c,x_c-r_ins,y_c,180,3)
            mi_clearselected;
            % specify coating material
            mi_addblocklabel(x_c-(r_ins+r_ext)/2,y_c);
            mi_selectlabel(x_c-(r_ins+r_ext)/2,y_c);
            mi_setblockprop(sprintf('cond%d_insu',j),0,0.001,0,0,0,1)
            mi_clearselected
        end
    end
    
    % draw domain limits and apply boundary conditions
    mi_addboundprop('A0', 0, 0, 0, 0, 0, 0, 0, 0, 0)
    
    B =[...
        L -L L L; ... %vertical line from y=-L to y=L, x=L
        -L -L -L L; ... %vertical line from y=-L to y=L, x=-L
        -L L L L; ... %horizontal line from x=-L to x=L, y=L
        -L -L L -L; ... %horizontal line from x=-L to x=L, y=-L
        ];
    for b=1:size(B,1)
        xs=B(b,1);ys=B(b,2);
        xe=B(b,3);ye=B(b,4);
        mi_drawline(xs,ys,xe,ye);
        mi_selectsegment((xs+xe)/2,(ys+ye)/2);
        mi_setsegmentprop('A0',0,1,0,0);
    end
    
    num_interfaces=length(t);
    for b=1:num_interfaces-1
        mi_drawline(-L,t(b),L,t(b));
    end
    
    window_size=(max(Geom(:,2))-min(Geom(:,2)))*.8;
    x_center=mean(Geom(:,2));
    y_center=mean(Geom(:,3));
    mi_zoom(x_center-window_size,y_center-window_size,x_center+window_size,y_center+window_size);
    
    %% Run solver & post-process
    mi_saveas(fname);
    mi_analyze(1)
    mi_loadsolution;
    mo_zoom(x_center-window_size,y_center-window_size,x_center+window_size,y_center+window_size);
    mo_shownames();
    
    % Inspect conductors
%         for j=8
%             x_c = Geom(j,2);
%             y_c = Geom(j,3);
%             r_ext = Geom(j,5);
%             r_ins = Geom(j,8);
%             DELTA=1e-6;
%             r_contour=sum([r_ext r_ins],'omitnan')+DELTA;
%             mo_seteditmode('contour');
%             mo_clearcontour();
%             window_size=3*r_contour;
%             mo_zoom(x_c-window_size,y_c-window_size,x_c+window_size,y_c+window_size);
%             theta=[0:(pi/200):2*pi];
%             for jj=1:length(theta)
%                 x=x_c+r_contour*cos(theta(jj));
%                 y=y_c+r_contour*sin(theta(jj));
%                 mo_addcontour(x,y);
%                 pause(0.01);
%             end
%         end
    
    % Extract impedances
    for j=i:con
        x_c = Geom(j,2);
        y_c = Geom(j,3);
        r_in = Geom(j,4);
        r_ext = Geom(j,5);
        %sig_c=1/Geom(j,6)/1e6;
        %r_ins = Geom(j,8);
        %r_total=sum([r_ext r_ins],'omitnan');
        dx=(r_ext+r_in)/2;
        mo_selectblock(x_c-dx,y_c);
        
        % Type    Definition
        % 0       A · J
        % 1       A
        % 2       Magnetic field energy
        % 3       Hysteresis and/or lamination losses
        % 4       Resistive losses
        % 5       Block cross-section area
        % 6       Total losses
        % 7       Total current
        % 8       Integral of Bx (or Br) over block
        % 9       Integral of By (or Bz) over block
        % 10      Block volume
        % 11      x (or r) part of steady-state Lorentz force
        % 12      y (or z) part of steady-state Lorentz force
        % 13      x (or r) part of 2× Lorentz force
        % 14      y (or z) part of 2× Lorentz force
        % 15      Steady-state Lorentz torque
        % 16      2× component of Lorentz torque
        % 17      Magnetic field coenergy
        % 18      x (or r) part of steady-state weighted stress tensor force
        % 19      y (or z) part of steady-state weighted stress tensor force
        % 20      x (or r) part of 2× weighted stress tensor force
        % 21      y (or z) part of 2× weighted stress tensor force
        % 22      Steady-state weighted stress tensor torque
        % 23      2× component of weighted stress tensor torque
        % 24      R2 (i.e. moment of inertia / density)
        % 25      x (or r) part of 1× weighted stress tensor force
        % 26      y (or z) part of 1× weighted stress tensor force
        % 27      1× component of weighted stress tensor torque
        % 28      x (or r) part of 1× Lorentz force
        % 29      y (or z) part of 1× Lorentz force
        % 30      1× component of Lorentz torque
        
        % S = mo_blockintegral(5);
        % A = mo_blockintegral(1)/S; % this is the INTEGRAL across the surface S, so it must be removed
        %
        % I = mo_blockintegral(7);
        % J = 1i*w*sig_c*A;
        % Z = J/sig_c;

        maxRetries = 5; % Maximum number of retries
        retryDelay = 0.25; % Delay between retries in seconds
        retryCount = 0; % Initialize retry count

        while retryCount < maxRetries
            try
                if i == j
                    vals = mo_getcircuitproperties(sprintf('cond%d_source',j));
                else
                    vals = mo_getcircuitproperties(sprintf('cond%d_target',j));
                end

                if numel(vals) > 1
                    break; % Exit the loop if vals has more than 1 element
                else
                    error('vals has only one element'); % Trigger an error to enter the catch block
                end
            catch ME
                % Increment retry count
                retryCount = retryCount + 1;

                % If max retries reached, rethrow the error
                if retryCount >= maxRetries
                    error('Failed to get valid circuit properties after %d retries: %s', maxRetries, ME.message);
                end

                % Pause before retrying
                pause(retryDelay);
            end
        end

        % while true
        %     if i==j
        %         vals = mo_getcircuitproperties(sprintf('cond%d_source',j));
        %     else
        %         vals = mo_getcircuitproperties(sprintf('cond%d_target',j));
        %     end
        %     if numel(vals)>1
        %         break
        %     end
        % end
        Z(i,j)= vals(2); % 2 is the voltage drop which is equal to J/sigma
        if i~=j;Z(j,i)=Z(i,j);end %symmetry is beautiful
    end
    mi_close
    
end %of main loop

end % of function