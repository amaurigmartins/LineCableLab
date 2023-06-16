function [] = set_plot_params(cfg)
% cfg must be a vector with [sk col k_scaling w_corr h_corr k_width_height fnt_scaling]
if nargin < 1
    cfg = 'default';
end

if strcmp(cfg,'default')
    sk = 1; %this will multiply everything for screen display purposes
    col = 1; % 1 column
    w_corr = 1;
    h_corr = 1;  % height correction to avoid bizarre windows behavior
    k_scaling = 2;          % scaling factor of the figure
    k_width_height = 1;      % width:height ratio of the figure
    fnt_scaling = 1;  % scaling factor for fonts
else
    sk = cfg(1);
    col = cfg(2);
    k_scaling = cfg(3);
    w_corr = cfg(4);
    h_corr = cfg(5);
    k_width_height = cfg(6);
    fnt_scaling = cfg(7);
end

% IEEE Standard Figure Configuration - Version 1.0
% run this code before the plot command

%
% According to the standard of IEEE Transactions and Journals: 
% Times New Roman is the suggested font in labels. 
% For a singlepart figure, labels should be in 8 to 10 points,
% whereas for a multipart figure, labels should be in 8 points.
% Width: column width: 8.8 cm; page width: 18.1 cm.

%
% width & height of the figure
% (You need to plot a figure which has a width of (8.8 * k_scaling)
% in MATLAB, so that when you paste it into your paper, the width will be
% scalled down to 8.8 cm  which can guarantee a preferred clearness.

% before we start just scale everything again to fit the screen
k_scaling = k_scaling * sk;
k_width_height = k_width_height * sk;
fnt_scaling = fnt_scaling * sk;

% now make everything look good
width = w_corr * (col * 8.8 * k_scaling);
height = h_corr * ((col * 8.8 * k_scaling) / k_width_height);

%% figure margins
top = 0.5 * h_corr;  % normalized top margin
bottom = 1.5 * h_corr;	% normalized bottom margin
left = 1.5 * w_corr;	% normalized left margin
right = 1 * w_corr;  % normalized right margin

%% set default figure configurations
reset(0);
fontface='Times New Roman';
set(0,'defaultFigureUnits','centimeters');
set(0,'defaultFigurePosition',[0 0 width height]);
set(0,'defaultLineLineWidth',1*k_scaling);
set(0,'defaultAxesLineWidth',0.25*k_scaling);
set(0,'defaultAxesGridLineStyle',':');
set(0,'defaultAxesYGrid','on');
set(0,'defaultAxesXGrid','on');
set(0,'defaultAxesFontName',fontface);
set(0,'defaultAxesFontSize',12*fnt_scaling);
set(0,'defaultTextFontName',fontface);
set(0,'defaultTextFontSize',12*fnt_scaling);
set(0,'defaultLegendFontName',fontface);
set(0,'defaultLegendFontSize',12*fnt_scaling);
set(0,'defaultAxesUnits','normalized');
set(0,'defaultAxesPosition',[left/width bottom/height (width-left-right)/width  (height-bottom-top)/height]);
% set(0,'defaultAxesColorOrder',[0 0 0]);
set(0,'DefaultAxesColorOrder','remove')
set(0,'DefaultAxesColorOrder',get(0,'FactoryAxesColorOrder'))
set(0,'defaultAxesLineStyleOrder','remove')
set(0,'defaultAxesLineStyleOrder',get(0,'FactoryAxesLineStyleOrder'))
set(0,'defaultAxesLineStyleOrder','-')
set(0,'defaultAxesTickDir','out');
set(0,'defaultFigurePaperPositionMode','auto');
set(0,'defaultLegendLocation','best');
set(0,'defaultLegendBox','on');
set(0,'defaultLegendOrientation','vertical');


end