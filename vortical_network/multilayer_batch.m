close all; clc; clear;
% ------------------------------------------------------------------

% Load in data set
%load('~/IBPM/mat_files/gust_0.3125_50_20_4.mat','xn','yn','omg','x_body','y_body', 'fxb', 'fyb', 'vxb', 'vyb');

load('./gust_0.3125_50_25_4.mat','xn','yn','omg','x_body','y_body', 'fxb', 'fyb', 'vxb', 'vyb');

% mat_name = sprintf('AoA35_gust_0.3125.mat');
% 
% fprintf('Loading dataset: %s\n', mat_name)
% 
% 
% load_path = sprintf('~/APS/mat/%s', mat_name);
% 
% load(load_path,'xn','yn','omg','x_body','y_body', 'fxb', 'fyb', 'vxb', 'vyb');

figure(10);
plot(x_body(end,:))


% ------------------------------------------------------------------
load colormap_latest.mat;
dx = xn(1,2) - xn(1,1);
dy = yn(2,1) - yn(1,1);
nt = 166; % number of times to loop over
xc_store = zeros(8,nt);
yc_store = zeros(8,nt);
Gammac_store = zeros(8,nt);
cmax_w = 4;
clev_w = 48; % use even #
clevs = linspace( -cmax_w, cmax_w, clev_w );
range = [-0.5 8.7 -2.5 2.5];
% ------------------------------------------------------------------

global start stop step dt it plot_dim top_old

start = 500;
stop = 501; %50
step = 1;
dt = 0.0025*20;

n_step = stop-start + 1;

top_old = 0;


% % % figure(99); hold all, axis off
% % % plot_dim = round(sqrt((stop - start)/step+1));
% % % %tiledlayout(8, 8,'TileSpacing','Compact');
% % % title(['time = ' num2str( it*dt )])

U_ref = 1;
N = size(x_body(:,1),1);
f = 1;
E = 1;
I = 1;
c = 1;

delta = sqrt(0.1*(c*(5/N)*U_ref/(f*E*I*10^6))^2);


% % %% Define what we want to save
% % 
% % wb = zeros(66, n_step);
% % A_fluid = zeros(10,10, n_step);
% % A_struct = zeros(10,10, n_step);
% % A_supra = zeros(20,20, n_step);
% % 
% % itr = 1;


for it = start:step:stop % must start from at least iteration 2!!!
    
    fprintf('\nIteration: %i\n', it)
    
    top = 0.01;
    
    vort = omg(:,:,it);
    xb = x_body(:,it);
    yb = y_body(:,it);
    
    % ---------------------------------------------------------------
    % ---------------------- Fluid Vorticity ------------------------
    % ---------------------------------------------------------------

    fprintf('! ---- Fluid Elements ---- !\n')

    % Define hyper-parameters
    gamma = 0.75;

    vort = omg(:,:,it);
    xb = x_body(:,it);
    yb = y_body(:,it);

    % fix large abs values
    wn = vort;
    wn( wn > cmax_w ) = cmax_w;
    wn( wn < -cmax_w ) = -cmax_w;

    % truncate the values of gamma and corrosponding x and y-position
    Gamma_v  = vort(abs(vort(:)) > 0.01*max(abs(vort(:))))*dx*dy;
    xfil    = xn(abs(vort(:)) > 0.01*max(abs(vort(:))));
    yfil    = yn(abs(vort(:)) > 0.01*max(abs(vort(:))));

    [Ci_fluid,A_fluid,xc_fluid,yc_fluid,Gamma_c_fluid] = plots_n_things(Gamma_v, xfil, yfil, gamma, 3, xb, yb, vort, xn, yn, 0);
    
    
    eval(['ADJ_fluid_', int2str(it),' = A_fluid']);
    % ---------------------- PLOTTING -----------------------------
% % % %     fprintf('Plotting fluid communities \n');
% % % % 
% % % %     figure(90), clf, axis off
% % % %     
% % % %     plt_dims = [-0.5 10.75 -3 2.5];
% % % %     c_dims = [-1 0];
% % % %     
% % % %     % Plot countours of vorticity
% % % %     ax1 = axes; hold on;
% % % %     contour(xn,yn,vort,linspace(-20,20),'linewidth',1.2)
% % % %     axis off, axis equal, axis(plt_dims), 
% % % %     colormap(ax1,parula);%colorbar;
% % % % 
% % % %     % Plot vortical community regions
% % % %     ax2 = axes; hold on;
% % % %     scatter(xfil,yfil,10,Ci_fluid,'filled','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8)
% % % %     axis off, axis equal, axis(plt_dims),
% % % %     colormap(ax2,parula);%colorbar;
% % % %     hold on;
% % % % 
% % % %     % Plot beam
% % % %     ax3 = axes; hold on;
% % % %     plot(xb,yb,'linewidth',5, 'color', 'b')
% % % %     axis off, axis equal, axis(plt_dims), 
% % % % 
% % % %     
% % % %     % Plot vortical community structure
% % % %     ax4 = axes; hold on;
% % % %     Node_weight = plotgraph(xc_fluid,yc_fluid,A_fluid); % plotgraph
% % % %     h = hot(512);
% % % %     colormap(ax4,flipud(h(1:end,:)));%colorbar;
% % % %     axis off, axis equal, axis(plt_dims)
% % % %     %caxis(ax3, c_dims);
% % % % 
% % % %     % Plot nodes
% % % %     for j = 1:length(xc_fluid)
% % % %         plot(xc_fluid(j), yc_fluid(j),'o','Color','k','MarkerFaceColor','k','MarkerSize',30*abs(log10(Node_weight(j))),'linewidth',1.7);
% % % %     end
% % % % 
% % % %     axis off, axis equal, axis(plt_dims);
% % % %     hold off
% % % %     
% % % %     % Plot adjacency matrix
% % % %     plot_adj(A_fluid, 56);
    
    % ---------------------------------------------------------------
    % -------------------- Structure Vortexlets ---------------------
    % ---------------------------------------------------------------

    fprintf('! ---- Structure Elements ---- !\n')
    
    % Define parameters
    delta = 0.001;
    gamma = 0.75;
    alpha = 50; % angle of attack in degrees
    fig_start = 1;
    
    % Solve for vortexlets on the structure
    [wb] = bound_vortexlet2(Gamma_c_fluid, xc_fluid, yc_fluid, vxb, vyb, xb, yb, it, delta, alpha); 
      %Mb, Mw, vel_control, pos_control, n,                                                           %(wn          , xn      , yn      , vxb, vyb, xb, yb, it, delta, alpha)
    
    eval(['wb_', int2str(it),' = wb']);
    
    
    Gamma_v = wb;
    xfil    = xb;
    yfil    = yb;

    [Ci,A_struct,xc_struct,yc_struct,Gamma_c_struct] = plots_n_things(wb, xb, yb, gamma, 4, xb, yb, vort, xn, yn, 1);
    
    eval(['ADJ_struct_', int2str(it),' = A_struct']);
    
    
% % % %     % ---------------------- PLOTTING -----------------------------
% % % %     fprintf('Plotting structural communities \n');
% % % %     
% % % %     plt_dims = [-0.25 1 -0.8750 0.3750];
% % % %     c_dims = [-1 0];
% % % %     
% % % %     % Convert color code to 1-by-3 RGB array (0~1 each)
% % % %     str = '#707070';
% % % %     color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
% % % % 
% % % %     figure(13), clf, axis off
% % % %        
% % % %     % Plot vortical community regions
% % % %     ax0 = axes; hold on;
% % % %     plot(xb,yb,'Color', color,'linewidth',5);
% % % %     axis off, axis equal, axis(plt_dims), 
% % % %     colormap(ax0,parula);%colorbar; 
% % % %     scatter(xfil,yfil, abs(Gamma_v)*4000,Ci,'filled')
% % % %     axis off, axis equal, axis(plt_dims) 
% % % %     
% % % %     ax2 = axes; hold on;
% % % %     Node_weight = plotgraph(xc_struct,yc_struct,A_struct); % plotgraph
% % % %     h = hot(512);
% % % %     colormap(ax2,flipud(h(1:end,:)));%colorbar;
% % % %     %caxis(ax2,c_dims);
% % % %     
% % % %     for j = 1:length(xc_struct)
% % % %         plot(xc_struct(j), yc_struct(j),'o','Color','k','MarkerFaceColor','k','MarkerSize',50*abs(log10(Node_weight(j))),'linewidth',1.7);
% % % %     end
% % % %     
% % % %     axis off, axis equal, axis(plt_dims);
% % % % 
% % % %     hold off
% % % %     
% % % %     % Plot adjacency matrix
% % % %     plot_adj(A_struct, 65);
    
    % ---------------------------------------------------------------
    % ---------------------- Supra-Adjacency ------------------------
    % ---------------------------------------------------------------

    A = zeros(length(A_struct)+length(A_fluid));
    pts = length(A);
    
    Gamma = cat(1, Gamma_c_struct, Gamma_c_fluid);
    Gamma = Gamma/(2*pi);
        
    xc = cat(1, xc_struct, xc_fluid);
    yc = cat(1, yc_struct, yc_fluid);
    r(1:pts, 1:2) = [xc yc]; % position vector matrix, r = (xi,yj)

    % Biot-savart law
    for i = 1:pts  
        dist_ij = r(1:pts,:) - ones(pts,1)*r(i,:);  % distance vector, r_{i->j}
        dist_mag = sqrt(sum(dist_ij.^2,2));         % ||r_{i->j}||
        dist_ij = 1./dist_mag;                      % distance ratio
        dist_ij(isinf(dist_ij)) = 0;                % update r_{i->i} = 0
        u_ij = abs(Gamma(i)) * dist_ij;             % u_{i->j}
        u_ji = abs(Gamma) .* dist_ij;               % u_{j->i}
        A(i,:) = u_ij;
    end
    
    eval(['ADJ_supra_', int2str(it),' = A']);
    % ---------------------- PLOTTING -----------------------------
    fprintf('Plotting communities \n');
    % --------------------- Reload Fluid ----------------------------
    % fix large abs values
    wn = vort;
    wn( wn > cmax_w ) = cmax_w;
    wn( wn < -cmax_w ) = -cmax_w;

    % truncate the values of gamma and corrosponding x and y-position
    Gamma_v  = vort(abs(vort(:)) > 0.01*max(abs(vort(:))))*dx*dy;
    xfil    = xn(abs(vort(:)) > 0.01*max(abs(vort(:))));
    yfil    = yn(abs(vort(:)) > 0.01*max(abs(vort(:))));
    
    % --------------------- plots start ----------------------------

    figure(66), clf, axis off
    
    plt_dims = [-0.5 10.75 -3 2.5];
    c_dims = [-1 0];
    
    % Plot countours of vorticity
    ax1 = axes; hold on;
    contour(xn,yn,vort,linspace(-20,20),'linewidth',1.2)
    axis off, axis equal, axis(plt_dims), 
    colormap(ax1,parula);%colorbar;

    % Plot vortical community regions
    ax2 = axes; hold on;
    scatter(xfil,yfil,10,Ci_fluid,'filled','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8)
    axis off, axis equal, axis(plt_dims),
    colormap(ax2,parula);%colorbar;
    hold on;

    % Convert color code to 1-by-3 RGB array (0~1 each)
    str = '#707070';
    color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

    % Plot vortical community regions
    %plot(xb,yb,'Color', color,'linewidth',5);
    ax3 = axes; hold on;
    Node_weight = plotgraph(xc,yc,A); % plotgraph
    h = hot(512);
    colormap(ax3,flipud(h(1:end,:)));%colorbar;
    axis off, axis equal, axis(plt_dims),
    %caxis(ax3,c_dims);

    % Plot nodes
    for j = 1:length(xc)
        plot(xc(j), yc(j),'o','Color','k','MarkerFaceColor','k','MarkerSize',30*abs(log10(Node_weight(j))),'linewidth',1.7);
    end

    axis off, axis equal, axis(plt_dims);
    hold off

    % Plot supra-adjacency matrix
    plot_adj(A, 9);
end

function [Ci,A_g,xc,yc,Gamma_c] = plots_n_things(Gamma_v, xfil, yfil, gamma, fig_start, xb, yb, vort, xn, yn, plot_or_not)
    % vortical communities                      (Gamma_v, xfil, yfil, gamma, 3,         xb, yb, vort, xn, yn, 0)
    fprintf('Running community detection \n');  

    % run vortical community detection
    [Ci,xc,yc,Gamma_c,~] = vortical_community(xfil,yfil,Gamma_v,1, gamma, 1); % vortical community
                                            
    [aa,bb] = sort(xc);
    xc = xc(bb);
    yc = yc(bb);
    Gamma_c = Gamma_c(bb);

    A_g = adjacency_mat(xc,yc,Gamma_c,1,'alg'); % adjacency matrix
end


function [] = plot_adj(ADJ, fig_start)

    global start stop step dt it top_old

    % Plot adjacency matrix
    %title(['time = ' num2str( it*dt )])

    %subplot(plot_dim,plot_dim,it/step,'Padding','Compact')

    top = max(max(ADJ));

    if top > top_old
        top_old = top;
    end

    %nexttile
    figure(fig_start);
    imagesc(ADJ);
    h = hot(512);
    colormap(gca,flipud(h(1:end,:)));
    caxis([0 1]);
    axis equal ij tight, axis off;

% %     if it == stop
% %         colorbar
% %     end
    %colorbar
    %title(['time = ' num2str( it*dt )])

% % %     N = length(ADJ);
% % %     x = repmat(1:N,N,1); % generate x-coordinates
% % %     y = x'; % generate y-coordinates
% % %     % Generate Labels
% % %     t = num2cell(ADJ); % extact values into cells
% % %     t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string
% % %     % Draw Image and Label Pixels
% % %     %%%%imagesc(M)
% % %     text(x(:), y(:), t, 'HorizontalAlignment', 'Center')


    %exportgraphics(gcf,['adj',num2str(it),'.png'],'Resolution',500);
end