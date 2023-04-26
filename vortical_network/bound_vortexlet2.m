function [wb] = bound_vortexlet2(wn, xn, yn, vxb, vyb, xb, yb, it, delta, alpha)

    % delta = 0.25; % smoothing parameter (published 0.0025)
    u = [1 0]; % freestream velocity
    nc = size(vxb,1); % number of control points
    
    % -------------------- Normal Unit Vectors --------------------------
    n = zeros(nc, 2);

    n(1, 1) = cosd(90 - alpha); 
    n(1, 2) = sind(90 - alpha); 

    for ii = 2:nc

       alpha = -atand(yb(ii)/xb(ii));

       n(ii, 1) = cosd(90 - alpha); 
       n(ii, 2) = sind(90 - alpha); 
       
    end
    
    
    % Check normal vector
    
% % %     figure(111);
% % %     plot(xb,yb); hold on;
% % %     quiver(xb,yb,n(:,1),n(:,2))
% % %     fprintf('Normals check: %f\n', sum((n(:,1).^2 + n(:,2).^2).^(0.5))/nc)
    
    % --------------------------------------------------------------------
    u_plate = [vyb vxb]; % set the normal velocity at each vortexlet equal to the normal velocity of the fluid
% % % %     u = repmat(u,nc,1); % put U(u, v) in a vector* nc long
% % % %     vel_control = zeros(nc,1); %initilize vector nc+1 long
% % % %    
     % ----------------------------- V*n ----------------------------------
    for i = 1:nc
        %vel_control(i,1) = dot(u_plate(i)-u(i),n(i)); % set the normal velocity at each vortexlet equal to the normal velocity of the fluid
        %vel_control(i,1) = dot((u(i)),n(i)); % set the normal velocity at each vortexlet equal to the normal velocity of the fluid
        %vel_control(i,1) = dot((u_plate(i)),n(i)); % set the normal velocity at each vortexlet equal to the normal velocity of the fluid
        vel_control(i,1) = (u_plate(i,1)^2 + u_plate(i,2)^2)^(1/2);
    end
% % % % 

    % co-locate control points with every vortexlet on the wing surface
    pos_control = [xb(:), yb(:)];    
    
    % ------------------------------------------------------------------
    % --------------------- "Matrix Maker" -----------------------------
    % ------------------------------------------------------------------
    
    % Solve for bound vortex matrix and evaluate strengths
   
    pos_vortex = pos_control;
    Mb = ones(nc); 
    
    for j = 1:nc
        for ii = 1:nc%length(pos_control(:,1))
            disp = pos_control(j,:) - pos_vortex(ii,:); % displacement from vortexlet to control point
            rsqr = disp(1).^2 + disp(2).^2; % radius squared of the displacement
            Mb(j,ii) = sum(conj([-disp(2)/(2*pi*(rsqr+delta^2)), disp(1)/(2*pi*(rsqr+delta^2))]).*n(j,:));
        end 
    end
    
    % -------------------------- Second Matrix ---------------------------
    
    pos_vortex = [xn, yn];
    Mw = ones(nc, length(wn)); 

    for j = 1:nc
        for ii = 1:length(wn)
            disp = pos_control(j,:) - pos_vortex(ii,:); % displacement from vortexlet to control point
            rsqr = disp(1).^2 + disp(2).^2; % radius squared of the displacement
            Mw(j,ii) = sum(conj([-disp(2)/(2*pi*(rsqr+delta^2)),disp(1)/(2*pi*(rsqr+delta^2))]).*n(j,:));
        end 
    end     
    
    % --------------------------------------------------------------------
%     
% %     Mb(end,:) = 1;
% %     Mw(end,:) = 1;
% %     vel_control(end) = 0;

    
    % solve for the strength of vortexlets gamma = inv(M)*-V
    strength_bound = Mb\(vel_control.*n - Mw*wn); 
       
    circ_at_each_vortex(:,1) = zeros(nc,1); 
    
    for pos = 1:nc
        circ_at_each_vortex(pos,2) = sum(strength_bound(1:pos));
    end
    
    %wb = circ_at_each_vortex(:,2)/(2*pi); % scalling for Gamma in 2D flow
    wb = circ_at_each_vortex(:,2);
end