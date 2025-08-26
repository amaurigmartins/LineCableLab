function [ResultTable, LCC_matrix] = calcCouplingRegionsParams(src,tgt,span_length)
%% DESCRIÇÃO:
%       Função que faz a subdivisão e encontra os parâmetros das regiões de acoplamento entre circuito fonte e alvo.
%
% INPUT:
%       src [m]: matriz com as coordenadas [X Y] dos pontos do circuito fonte.
%       tgt [m]: matriz com as coordenadas [X Y] dos pontos do circuito alvo.
%       span_length [m]: constante para definir a cada quantos metros o circuito alvo será subdividido
%
% OUTPUT:
%       ResultTable = [x_src y_src x_int_tgt y_int_tgt d1 d2 Leq_src Leq_tgt theta idx_acop Dist_acc_tgt]
%           x_src [m]: coordenada X do ponto no circuito fonte;
%           y_src [m]: coordenada Y do ponto no circuito fonte;
%           x_int_tgt[m]: coordenada X do ponto de acoplamento mais próximo no circuito alvo;
%           y_int_tgt[m]: coordenada Y do ponto de acoplamento mais próximo no circuito alvo;
%           d1 [m]: distância 1 do ponto do circuito fonte ao ponto de acoplamentoo do circuito alvo;
%           d2 [m]: distância 2 do ponto do circuito fonte ao ponto de acoplamentoo do circuito alvo;
%           Leq_src [m]: Comprimento equivalente com referência no circuito fonte
%           Leq_tgt [m]: Comprimento equivalente com referência no circuito alvo
%           idx_acop [adimensional]: índice que representa o tipo de acoplamento existente, de acordo com:
%                   0 : acoplamento normal entre o vão do circuito fonte e
%                   circuito alvo
%                   1 : acoplamento subdividido dentro de um mesmo vão
%           Theta [rad]: ângulo entre os trechos do circuito fonte e alvo
%           Dist_acc_tgt[m]: Distância acumulada com relação ao circuito alvo
%
%       LCC_matrix = [deq Leq Dist_acc_tgt idx_acop Leq_tgt]
%           deq [m]: distância equivalente entre circuito fonte e alvo
%           Leq [m]: comprimento equivalente entre circuito fonte e alvo
%           Dist_acc_tgt[m]: Distância acumulada com relação ao circuito alvo
%           idx_acop [adimensional]: índice que representa o tipo de acoplamento existente, de acordo com:
%                   0 : acoplamento normal entre o vão do circuito fonte e
%                   circuito alvo
%                   1 : acoplamento subdividido dentro de um mesmo vão
%           Leq_tgt [m]: Comprimento equivalente com referência no circuito alvo

%% PLOT

% set_plot_params;
close all
debugthis_closer_points = true; % FIGURE 1 - MACRO-REGIÕES
debugthis = false; % FIGURE 2 = MICRO-REGIÕES

%% PRÉ-PROCESSAMENTO:

[src, tgt] = pre_processing_polyline(src,tgt);

% Subdivide o duto a cada span_length metros:
[tgt] = sub_line(src,tgt,span_length);

% Encontra onde o DT e a LT se cruzam
[x_int,y_int,idx_int] = polyxpoly(src(:,1),src(:,2),tgt(:,1),tgt(:,2));
P_int = [x_int y_int];

%% DETERMINANDO MACRO-REGIÕES DE ACOPLAMENTO:

macro_regions_precision = 0.1; % Define a precisão para subdividir o duto e encontrar os pontos da macro-região

[tgt_sub] = sub_line(src,tgt,macro_regions_precision);

[~, tgt_closer_points] = closer_points2(src,tgt_sub,P_int,debugthis_closer_points);

[minor_d, ~] = closer_points2(src,tgt,P_int,0);


%% DETERMINANDO OS PARÂMETROS DE ACOPLAMENTO:

% [coord_XY_f,D1D2_f,Leq_src_f,Leq_tgt_f,Theta_f,IDX_f,dist_acc_f] = calculate_sub_params(src,tgt,minor_d,tgt_closer_points,debugthis);
[coord_XY_f,D1D2_f,Leq_src_f,Leq_tgt_f,Theta_f,IDX_f,dist_acc_f,dist_acc_src_f] = calculate_sub_params(src,tgt,minor_d,tgt_closer_points,debugthis);

ResultTable = [coord_XY_f D1D2_f Leq_src_f Leq_tgt_f Theta_f IDX_f dist_acc_f dist_acc_src_f];

LCC_matrix = define_LCC_matrix(ResultTable,idx_int);


end


function [new_src, new_tgt] = pre_processing_polyline(src,tgt)


% Determine domain boundaries
[minx,maxx,miny,maxy]=dombound(src,tgt);

% Retira pontos duplicados:
new_src = unique(src,'rows','stable');
new_tgt = unique(tgt,'rows','stable');

% Verifica se a polyline ta invertida e resolve o problema:
[~,~,alldistances] = minor_distance([maxx 0],new_src);

if alldistances(1) < alldistances(end)
    new_src = flip(new_src,1);
else
    new_src = src;
end

[~,~,alldistances] = minor_distance([maxx 0],new_tgt);

if alldistances(1) < alldistances(end)
    new_tgt = flip(new_tgt,1);
else
    new_tgt = tgt;
end



end

function [ minx,maxx,miny,maxy ] = dombound( src, tgt )
% dombound: Find domain limits based on supplied polylines.

minx=min([src(:,1);tgt(:,1)]);
maxx=max([src(:,1);tgt(:,1)]);

miny=min([src(:,2);tgt(:,2)]);
maxy=max([src(:,2);tgt(:,2)]);

end


function [m,idx,alldistances] = minor_distance(P,line)

X = P(1,1);
Y = P(1,2);
alldistances= sqrt((X-line(:,1)).^2 + (Y-line(:,2)).^2);
[m,idx]=min(alldistances);

end

function [new_line] = sub_line(src,line,n_sub)
% Divide a linha a cada n_sub metros

src_row = size(src,1);

debugthis = true;

[dist_total,M_dist_acc] = dist_acum(line);

lim = dist_total/src_row;

if n_sub > lim
    n_sub = lim;
    dist = linspace(0,dist_total,round(dist_total/n_sub));
else
    dist = linspace(0,dist_total,(round(dist_total/n_sub)+1));
end

Y = line(:,2);
X = line(:,1);

Y1 = interp1(M_dist_acc,Y,dist);
X1 = interp1(M_dist_acc,X,dist);

X1 = transpose(X1);
Y1 = transpose(Y1);

new_line = [X1 Y1];

if debugthis
    figure(1)
    plot(X,Y,'b^',X1,Y1,'r*');
    legend('line','new line');
end


end


function [dist_acc,M_dist_acc,M_Leq] = dist_acum(line)
% Função que calcula o valor total de comprimento da linha

dist_acc = 0;
M_dist_acc(1) = 0;
M_Leq(1) = 0;
for i = 1:size(line,1)-1
    X1 = line(i,1);
    Y1 = line(i,2);
    X2 = line(i+1,1);
    Y2 = line(i+1,2);

    dist = sqrt((X1 - X2)^2 + (Y1-Y2)^2);

    M_Leq(i+1) = dist;

    dist_acc = dist_acc + dist;

    M_dist_acc(i+1) = dist_acc;
end
M_dist_acc = transpose(M_dist_acc);

end

function [minor_d, tgt_closer_points] = closer_points2(new_src,new_tgt,P_int,debugthis)

if size(new_src,1) ~= size(new_tgt,1)

    tgt_comp = new_tgt;

    new_tgt(end,:) = [];

    for i = 1: size(new_src,1)
        X = new_src(i,1);
        Y = new_src(i,2);

        if i == 1
            tgt_closer_points(i,:) = new_tgt(1,:);
            minor_d(i,1) = 1;

        elseif i == size(new_src,1)
            tgt_closer_points(i,:) = tgt_comp(end,:);
            minor_d(i,1) = size(tgt_comp,1);

        else

            alldistances= sqrt((X-new_tgt(:,1)).^2 + (Y-new_tgt(:,2)).^2);
            [m,idx]=min(alldistances);
            tgt_closer_points(i,:) = new_tgt(idx,:);

            if new_tgt(idx,:) ~= P_int
                tgt_bkp=new_tgt;
                new_tgt(1:idx(1),:) = [];
                if isempty(new_tgt);new_tgt=tgt_bkp;end;
            end

            [m,idx,alldistances] = minor_distance(tgt_closer_points(i,:),tgt_comp);
            minor_d(i,1) = idx(1);


        end

        if debugthis
            figure(1)
            plot(tgt_comp(:,1),tgt_comp(:,2),'b',new_src(:,1),new_src(:,2),'r')
            xlabel('X [m]')
            ylabel('Y [m]')
            title('Coupling regions')
            axis equal
            hold on
            plot([new_src(i,1) tgt_closer_points(i,1)],[new_src(i,2) tgt_closer_points(i,2)],'k',new_src(:,1),new_src(:,2),'k^',tgt_closer_points(:,1),tgt_closer_points(:,2),'r^')
            legend('Target','Source','Coupling region projection')
        end

    end

    minor_d = sort(minor_d);
else
    minor_d = [1:1:size(new_src)];
    tgt_closer_points = new_tgt;
end


end

function [coord_XY_f,D1D2_f,Leq_src_f,Leq_tgt_f,Theta_f,IDX_f,dist_acc_f,dist_acc_src_f] = calculate_sub_params(src,tgt,minor_d,tgt_closer_points,debugthis)

D1D2_f = [];
Leq_src_f = [];
Leq_tgt_f = [];
Theta_f = [];
IDX_f = [];
dist_acc_tgt = 0;
dist_acc_src = 0;
dist_acc_f = [];
coord_XY_f = [];
dist_acc_src_f = [];

for i = 1: size(src,1) - 1
    src_span = src(i:i+1,:);
    if minor_d(i) == minor_d(i+1)
        tgt_span = [tgt(minor_d(i),:); tgt(minor_d(i),:)];
    else
        tgt_span = tgt(minor_d(i):minor_d(i+1),:);
    end

    N = size(tgt_span,1); % Número de micro regiões
    X_src = linspace(src_span(1,1),src_span(2,1),N);
    Y_src = linspace(src_span(1,2),src_span(2,2),N);

    X_src = transpose(X_src);
    Y_src = transpose(Y_src);

    src_span_sub = [X_src Y_src];

    if debugthis
        figure(2)
        plot(src_span(:,1),src_span(:,2),'r',tgt_span(:,1),tgt_span(:,2),'b');
        xlabel('X [m]')
        ylabel('Y [m]')
        hold on
    end

    for j = 1: N - 1
        D(j,1) = dist_2p(src_span_sub(j,:),tgt_span(j,:));
        Leq_src(j,1) = dist_2p(src_span_sub(j,:),src_span_sub(j+1,:));
        Leq_tgt(j,1) = dist_2p(tgt_span(j,:),tgt_span(j+1,:));
        Theta(j,1) = make_angle(src_span_sub(j:j+1,:),tgt_span(j:j+1,:));
        dist_acc_tgt = dist_acc_tgt + Leq_tgt(j,1);
        dist_acc_f = [dist_acc_f;dist_acc_tgt];
        dist_acc_src = dist_acc_src + Leq_src(j,1);
        dist_acc_src_f = [dist_acc_src_f;dist_acc_src];

        if debugthis
            figure(2)
            plot([src_span_sub(j,1); tgt_span(j,1)],[src_span_sub(j,2); tgt_span(j,2)],'k--',src_span(:,1),src_span(:,2),'k^',tgt_closer_points(:,1),tgt_closer_points(:,2),'r^',src_span_sub(:,1),src_span_sub(:,2),'rx',tgt_span(:,1),tgt_span(:,2),'bx');
            legend('Source','Target','Projection')
            title('Subdivisions')
            axis equal
            hold on
        end

    end

    D(N,1) = dist_2p(src_span_sub(end,:),tgt_span(end,:));
    D1D2 = [D(1:N-1) D(2:N)];

    if N ~= 1
        %         if i == size(src,1) - 1
        %             IDX = zeros(N-1,1);
        %             IDX(1) = 1;
        %             IDX(end) = 1;
        %         elseif i == 1 && size(src,1) ~= size(tgt,1)
        %             IDX = zeros(N-1,1);
        %         else
        IDX = zeros(N-1,1);
        IDX(end) = 1;
        %         end
    else
        IDX = [];
    end

    coord_XY = [src_span_sub(1:end-1,:) tgt_span(1:end-1,:)];

    % Montando a lista final com todos os parâmetros:
    D1D2_f = [D1D2_f;D1D2];
    Leq_src_f = [Leq_src_f;Leq_src];
    Leq_tgt_f = [Leq_tgt_f;Leq_tgt];
    Theta_f = [Theta_f;Theta];
    IDX_f = [IDX_f;IDX];
    coord_XY_f = [coord_XY_f;coord_XY];

    % Limpando as variáveis:
    Leq_src = [];
    Leq_tgt = [];
    Theta = [];
    D = [];
    IDX = [];

end


end


function D = dist_2p(P1,P2)
% Determina a distancia entre os pontos P1 e P2

D = norm(P1-P2);

end

function LCC_matrix = define_LCC_matrix(TR,idx_int)
% Função que gera a matrix de parâmetros para os blocos LCC
%
% INPUT:
%       TR = [x_src y_src x_int_tgt y_int_tgt d1 d2 Leq_src Leq_tgt theta idx_acop Dist_acc_tgt]
%           x_src [m]: coordenada X do ponto no circuito fonte;
%           y_src [m]: coordenada Y do ponto no circuito fonte;
%           x_int_tgt[m]: coordenada X do ponto de acoplamento mais próximo no circuito alvo;
%           y_int_tgt[m]: coordenada Y do ponto de acoplamento mais próximo no circuito alvo;
%           d1 [m]: distância 1 do ponto do circuito fonte ao ponto de acoplamentoo do circuito alvo;
%           d2 [m]: distância 2 do ponto do circuito fonte ao ponto de acoplamentoo do circuito alvo;
%           Leq_src [m]: Comprimento equivalente com referência no circuito fonte
%           Leq_tgt [m]: Comprimento equivalente com referência no circuito alvo
%           idx_acop [adimensional]: índice que representa o tipo de acoplamento existente, de acordo com:
%                   0 : acoplamento normal entre o vão do circuito fonte e
%                   circuito alvo
%                   1 : acoplamento subdividido dentro de um mesmo vão
%           Theta [rad]: ângulo entre os trechos do circuito fonte e alvo
%           Dist_acc_tgt[m]: Distância acumulada com relação ao circuito alvo
%       idx_int: linha em que ocorre o cruzamento entre o circuito fonte e alvo
%
%
% OUTPUT:
%       LCC_matrix = [deq Leq Dist_acc_tgt idx_acop Leq_tgt]
%           deq [m]: distância equivalente entre circuito fonte e alvo
%           Leq [m]: comprimento equivalente entre circuito fonte e alvo
%           Dist_acc_tgt[m]: Distância acumulada com relação ao circuito alvo
%           idx_acop [adimensional]: índice que representa o tipo de acoplamento existente, de acordo com:
%                   0 : acoplamento normal entre o vão do circuito fonte e
%                   circuito alvo
%                   1 : acoplamento subdividido dentro de um mesmo vão
%           Leq_tgt [m]: Comprimento equivalente com referência no circuito alvo


deq = (TR(:,5) + TR(:,6))./2;

Leq = sqrt(TR(:,7).*TR(:,8)).*abs(cos(TR(:,9)));

LCC_matrix = [deq Leq TR(:,11) TR(:,10) TR(:,8) TR(:,12)];

if ~isempty(idx_int)
    for i = 1: size(idx_int,1)
        LCC_matrix(idx_int(i,2)+1:end,1) = -LCC_matrix(idx_int(i,2)+1:end,1);
    end
end

end

function [Theta] = make_angle(src,tgt)
% Função que calcula o angulo(radianos) entre os trechos

if size(src,1) == 2
    v_src = [src(2,1:2) 0] - [src(1,1:2) 0];
    v_tgt = [tgt(2,1:2) 0] - [tgt(1,1:2) 0];
    Theta(1,1) = atan2(norm(cross(v_src, v_tgt)), dot(v_src, v_tgt));
else
    for i = 1 : size(src,1)-1

        v_src = [src(i+1,1:2) 0] - [src(i,1:2) 0];
        v_tgt = [tgt(i+1,1:2) 0] - [tgt(i,1:2) 0];
        Theta(i,1) = atan2(norm(cross(v_src, v_tgt)), dot(v_src, v_tgt));

    end
end

end