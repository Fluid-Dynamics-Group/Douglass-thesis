function Node_weight = plotgraph(A_x,A_y,A_g)
%%
% Network representation of flow field
% Date Created: 06/30/2016
% Date Modified: 01/05/2017, general flow field plot(plot_flow)
% By MGM
% Reference: M. Gopalakrishnan Meena, A. G. Nair & K. Taira
%    "Network Community-Based Model Reduction for Vortical Flows",
%    Physical Review E, 97, 063103, 2018 (doi.org/10.1103/PhysRevE.97.063103)
%%
% plotgraph_flow(itr,A_x,A_y,A_g,Xlim,Ylim)
% Function to visualize the network of a vortex system in a flow field
% representation
% body_data : (1)=body type(1=cylinder, 2=airfoil); (2)=itr; (3-6)=airfoil
%             body info(3=m, 4=L, 5=h, 6=t)
% A_x       : x coordinates of vertices
% A_y       : y coordinates of vertices
% A_g       : adjacency matrix
% Xlim      : x limits of the field
% Ylim      : y limits of the field

%% Load network details of current snapshot
x = A_x; y = A_y;
pts = length(x)-1;
n_cor = zeros(pts+1,2);     % coordinates of vortices and cylinder
for j = 1:pts+1
    n_cor(j,:) = [x(j),y(j)];
end
A_g_nor = A_g/max(max(A_g));
Node_weight = (sum(A_g_nor,1)' + sum(A_g_nor,2))*2;
%% Vortex centers
% for j = 1:pts+1
%     plot(ax1,n_cor(j,1), n_cor(j,2),'o','Color',color(2,:),...%'k',...
%         'MarkerFaceColor',color(2,:),...
%         'MarkerSize',30 * abs(log10(Node_weight(j))),'linewidth',1.7);      
% end
%ax23 = axes;
%% Flow field representation of network
for j = 1:pts+1
    if j/2 == 0 % change parabola orientation for odd and even starting points
        dir = -1;
    else
        dir = 1;
    end
    for k = j+1:pts+1
        edge_color = A_g(j,k);
        trans_color = A_g_nor(j,k);
        x_pts = [n_cor(j,1), n_cor(k,1)];
        y_pts = [n_cor(j,2), n_cor(k,2)];
        % create parabolic lines between points
        % y = ax^2 + bx + c: y(1) and y(2) - 2 eqns.
        % y' = 2ax + b = slope of parabola: y'(1) - 1 eqn.
        % take y'(1), y(1) and y(2) - solve 3 equations to obtain a,b,c
        A_p = [2*x_pts(1),    1,          0;
            x_pts(1)^2,    x_pts(1),   1;
            x_pts(2)^2,    x_pts(2),   1];
        B_p = [ dir*1.0; % slope of parabola changes sign alternatively
            y_pts(1);
            y_pts(2)];
        C_p = A_p\B_p;
        a = C_p(1); b = C_p(2); c = C_p(3); % coefficients of parabola equation
        yp = @(xp) a*xp.^2 + b*xp + c;
        x_co = linspace(x_pts(1),x_pts(2));
        y_co = yp(x_co);
        y_co(end) = NaN; % to remove line between first and last point
        patch(x_co, y_co, log10( edge_color*ones(1,length(x_co)) ), ...
            'EdgeColor', 'flat', 'FaceColor','none','LineWidth',4.5,...
            'EdgeAlpha',trans_color)
        dir = -dir;
    end

end
end     % end function

