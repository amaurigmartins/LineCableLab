function [Y]=Y_femm_slf_mut(Geom,soilFD,k,freq,con,basename)

Y = zeros(con,con);

try; closefemm; end

for i=1:con
    %% Open a FEMM instance
    fclose('all');
    
    try; closefemm; end
    fname=sprintf('%s_c%d.fec',basename,i);


    global HandleToFEMM;
    openfemm;
    hand1=HandleToFEMM;
    % create new document
    DOCTYPE=3; %0 for magnetics, 1 for electrostatics, 3 for current flow
    newdocument(DOCTYPE);
    % main_minimize;
    
    %% Set problem solver
    w=2*pi*freq;
    minangle=33.8;
    depth=1;
    ci_smartmesh(1)
    ci_probdef('meters','planar',freq,1e-12,depth,minangle);

    %ci_setgrid(1,'cart')
    eps0=8.8541878128e-12;
    mu0=4*pi*1e-7;
    
    %% Define energizations
    Vsrc=1;
    for j=1:con
        if j==i
            ci_addconductorprop(sprintf('cond%d_source',j), Vsrc, 0, 1)
        else
            ci_addconductorprop(sprintf('cond%d_target',j), 0, 0, 1)
        end
    end

    %% Define materials
    % Air and earth layers
    ci_addmaterial('air',0,0,1,1,0,0) % air layer
    if isfield(soilFD,'layer')
        % will handle this when I handle this
    else
        % m_g=soilFD.m_g/mu0;
        eps_g=soilFD.erg_total(k);
        sig_g=soilFD.sigma_g_total(k);
        ci_addmaterial('layer1', sig_g, sig_g, eps_g, eps_g, 0, 0)
    end
    
    % Conductors and insulation
    for j=1:con
        % m_c=Geom(j,7);
        sig_c=1/Geom(j,6);
        ci_addmaterial(sprintf('cond%d_core',j), sig_c, sig_c, 1, 1, 0, 0)
        if ~isnan(Geom(j,8))
            m_i=Geom(j,9);
            sig_i=0;
            eps_i=Geom(j,10);
            ci_addmaterial(sprintf('cond%d_insu',j), sig_i, sig_i, eps_i, eps_i, 0, 0) % conductor insulation
        end
    end
    
    %% Position objects in the 2D grid space
    L=20e3; % domain size
    
    % Air layer
    ci_addblocklabel(0,L/2);
    ci_selectlabel(0,L/2);
    ci_setblockprop('air',1,1,0);
    ci_clearselected;

    
    % Earth layers
    if isfield(soilFD,'layer')
        % will handle this when I handle this
    else
        t=[0 -L]; % soil layer interfaces
        ci_addblocklabel(0,(t(end)-t(1))/2);
        ci_selectlabel(0,(t(end)-t(1))/2);
        ci_setblockprop('layer1',1,1,0);
        ci_clearselected;
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
            ci_drawarc(x_c-r_in,y_c,x_c+r_in,y_c,180,3)
            ci_drawarc(x_c+r_in,y_c,x_c-r_in,y_c,180,3)
            ci_addblocklabel(x_c,y_c);
            ci_selectlabel(x_c,y_c);
            ci_setblockprop('air',1,0,0);
        end
        ci_clearselected

        % draw outer circle
        ci_drawarc(x_c-r_ext,y_c,x_c+r_ext,y_c,180,3)
        ci_drawarc(x_c+r_ext,y_c,x_c-r_ext,y_c,180,3)
        ci_addblocklabel(x_c+(r_in+r_ext)/2,y_c);
        ci_selectlabel(x_c+(r_in+r_ext)/2,y_c);
        if j==i
            ci_setblockprop(sprintf('cond%d_core',j),0.001,0,0);
            ci_clearselected;
            ci_selectarcsegment(x_c-(r_ext+r_in)/2,y_c-r_ext);
            ci_selectarcsegment(x_c-(r_ext+r_in)/2,y_c+r_ext);
            ci_setarcsegmentprop(.1, 0, 0, 0, sprintf('cond%d_source',j));
            %ci_setarcsegmentprop(maxsegdeg, "propname", hide, group, "inconductor") 
        else
            ci_setblockprop(sprintf('cond%d_core',j),0.001,0,0);
            ci_clearselected;
            %ci_seteditmode('arcsegments')
            ci_selectarcsegment(x_c-(r_ext+r_in)/2,y_c-r_ext);
            ci_selectarcsegment(x_c-(r_ext+r_in)/2,y_c+r_ext);
            ci_setarcsegmentprop(.1, 0, 0, 0, sprintf('cond%d_target',j));
        end
        
        % draw insulation
        if ~isnan(r_ins)
            ci_drawarc(x_c-r_ins,y_c,x_c+r_ins,y_c,180,3)
            ci_drawarc(x_c+r_ins,y_c,x_c-r_ins,y_c,180,3)
            ci_clearselected;
            % specify coating material
            ci_addblocklabel(x_c-(r_ins+r_ext)/2,y_c);
            ci_selectlabel(x_c-(r_ins+r_ext)/2,y_c);
            ci_setblockprop(sprintf('cond%d_insu',j),0.001,0,0);
            ci_clearselected
            ci_selectarcsegment(x_c-(r_ext+r_ins)/2,y_c-r_ins);
            ci_selectarcsegment(x_c-(r_ext+r_ins)/2,y_c+r_ins);
            ci_setarcsegmentprop(.1, 0, 0, 0, 0);
        end
    end
    
    % draw domain limits and apply boundary conditions
    ci_addboundprop('A0', 0, 0, 0, 0);
    
    B =[...
        L -L L L; ... %vertical line from y=-L to y=L, x=L
        -L -L -L L; ... %vertical line from y=-L to y=L, x=-L
        -L L L L; ... %horizontal line from x=-L to x=L, y=L
        -L -L L -L; ... %horizontal line from x=-L to x=L, y=-L
        ];
    for b=1:size(B,1)
        xs=B(b,1);ys=B(b,2);
        xe=B(b,3);ye=B(b,4);
        ci_drawline(xs,ys,xe,ye);
        ci_selectsegment((xs+xe)/2,(ys+ye)/2);
        ci_setsegmentprop('A0',0,1,0,0,0);
    end
    
    num_interfaces=length(t);
    for b=1:num_interfaces-1
        ci_drawline(-L,t(b),L,t(b));
    end
    
    window_size=(max(Geom(:,2))-min(Geom(:,2)))*.8;
    x_center=mean(Geom(:,2));
    y_center=mean(Geom(:,3));
    ci_zoom(x_center-window_size,y_center-window_size,x_center+window_size,y_center+window_size);
    
    %% Run solver & post-process
    ci_saveas(fname);
    ci_analyze(1)
    ci_loadsolution;
    co_zoom(x_center-window_size,y_center-window_size,x_center+window_size,y_center+window_size);
    co_shownames();
    
    % Inspect conductors
    %     for j=1
    %         x_c = Geom(j,2);
    %         y_c = Geom(j,3);
    %         r_ext = Geom(j,5);
    %         r_ins = Geom(j,8);
    %         DELTA=1e-3;
    %         r_contour=sum([r_ext r_ins],'omitnan')+DELTA;
    %         co_seteditmode('contour');
    %         co_clearcontour();
    %         window_size=3*r_contour;
    %         co_zoom(x_c-window_size,y_c-window_size,x_c+window_size,y_c+window_size);
    %         theta=[0:(pi/200):2*pi];
    %         for jj=1:length(theta)
    %             x=x_c+r_contour*cos(theta(jj));
    %             y=y_c+r_contour*sin(theta(jj));
    %             co_addcontour(x,y);
    %             pause(0.01);
    %         end
    %     end
    
    % Extract admittances
    for j=i:con
        x_c = Geom(j,2);
        y_c = Geom(j,3);
        r_in = Geom(j,4);
        r_ext = Geom(j,5);
        r_ins = Geom(j,8);
        % DELTA=1e-3;
        % r_contour=sum([r_ext r_ins],'omitnan')+DELTA;
        % co_seteditmode('contour');
        % co_clearcontour();
        % window_size=3*r_contour;
        % co_zoom(x_c-window_size,y_c-window_size,x_c+window_size,y_c+window_size);
        % theta=[0:(pi/200):2*pi];
        % for jj=1:length(theta)
        %     x=x_c+r_contour*cos(theta(jj));
        %     y=y_c+r_contour*sin(theta(jj));
        %     co_addcontour(x,y);
        %     pause(0.01);
        % end
        % Y=co_lineintegral(1);

        maxRetries = 5; % Maximum number of retries
        retryDelay = 0.25; % Delay between retries in seconds
        retryCount = 0; % Initialize retry count

        while retryCount < maxRetries
            try
                if i == j
                    vals = co_getconductorproperties(sprintf('cond%d_source', j));
                else
                    vals = co_getconductorproperties(sprintf('cond%d_target', j));
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
                    error('Failed to get valid conductor properties after %d retries: %s', maxRetries, ME.message);
                end

                % Pause before retrying
                pause(retryDelay);
            end
        end

        Y(i,j)= vals(2); 
        if i~=j;Y(j,i)=Y(i,j);end %symmetry is beautiful
    end
    ci_close
    
end %of main loop

end % of function