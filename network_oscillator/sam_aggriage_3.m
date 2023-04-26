clear;clc;close all;

% Load data
% ------------------------------------------------------------------
load data_full_sam.mat;   % zprime_norm_amp_1 zprime_norm_dot_amp_1 54900x4
inter = 20;               % time step between snapshots
dt    = 0.0025*inter;     % time interval
nt    = 25000/inter;      % number of snapshots
time  = 0:dt:dt*(nt-1);   % set time vector
nosc  = 7;
% load test_data_two.mat;  % zprime_test, zprimedot_test 117120x4

zrefd = (z_dot_test);

% training data - pert1
% ------------------------------------------------------------------
z_train     = z_amp_1;
z_dot_train = z_dot_amp_1;

Areg1           = regression_reconstruct(size(z_dot_train ,1), size(z_dot_train ,2), z_train, z_dot_train);
Lreg1           = Areg1 - diag(sum(Areg1,2));
z_reg_dot1      = zeros(size(z_test));

for i = 1:size(z_test,1)
    z_reg_dot1(i,:) = odefun(1, z_test(i,1:nosc), Areg1);
end

zrd1 = (z_reg_dot1);

% training data - pert2
% ------------------------------------------------------------------
z_train     = z_amp_2;
z_dot_train = z_dot_amp_2;

Areg2           = regression_reconstruct(size(z_dot_train ,1), size(z_dot_train ,2), z_train, z_dot_train);
Lreg2           = Areg2 - diag(sum(Areg2,2));
z_reg_dot2      = zeros(size(z_test));

for i = 1:size(z_test,1)
    z_reg_dot2(i,:) = odefun(1, z_test(i,1:nosc), Areg2);
end

zrd2 = (z_reg_dot2);

% % training data - pert3
% % ------------------------------------------------------------------
% zprime_train    = zprime_norm_amp_3;
% zprimedot_train = zprime_norm_dot_amp_3;
% 
% Areg3           = regression_reconstruct(size(zprimedot_train ,1), size(zprimedot_train ,2), zprime_train, zprimedot_train);
% Lreg3           = Areg3 - diag(sum(Areg3,2));
% z_reg_dot3      = zeros(size(zprime_test));
% 
% for i = 1:size(zprime_test,1)
%     z_reg_dot3(i,:) = odefun(1, zprime_test(i,1:nosc), Areg3);
% end
% 
% zrd3 = (z_reg_dot3);

% training data - combined
% ------------------------------------------------------------------
%load train_data_two.mat; % zprime_train, zprimedot_train 468480x4 (just all 4 pert cases vertically stacked, sepcifically the test data)
z_train     = z_train_all;
z_dot_train = z_dot_train_all;

bnds = 100000;

Areg      = regression_reconstruct(size(z_dot_train(1:bnds,:) ,1), size(z_dot_train(1:bnds,:) ,2), z_train(1:bnds,:), z_dot_train(1:bnds,:));
Lreg      = Areg - diag(sum(Areg,2));
z_reg_dot = zeros(size(z_test));

save("sam_aggrigate_L.mat", 'Lreg')

for i = 1:size(z_test,1)
    z_reg_dot(i,:) = odefun(1, z_test(i,1:nosc), Areg);
end

zrd = (z_reg_dot);


% ------------------------------------------------------------------
% Plotting
% ------------------------------------------------------------------


% Error plotting (Eqn 14)
% ------------------------------------------------------------------
r1 = sqrt( sum( abs((zrefd(:,:) - zrd1(:,:)).^2)) / numel(zrefd) )./(max(abs(zrefd(:,:))))*100;
r2 = sqrt( sum( abs((zrefd(:,:) - zrd2(:,:)).^2)) / numel(zrefd) )./(max(abs(zrefd(:,:))))*100;
% r3 = sqrt( sum( abs((zrefd(:,:) - zrd3(:,:)).^2)) / numel(zrefd) )./(max(abs(zrefd(:,:))))*100;
r  = sqrt( sum( abs((zrefd(:,:) - zrd(:,:)).^2 )) / numel(zrefd) )./(max(abs(zrefd(:,:))))*100;

figure;
subplot(221);hold all;
plot(r1,'o-','Color',color_save_new(1),'MarkerEdgeColor',color_save_new(1), 'MarkerFaceColor',color_save_new(1),'Linewidth',2);
plot(r2,'o-','Color',color_save_new(2),'MarkerEdgeColor',color_save_new(2), 'MarkerFaceColor',color_save_new(2),'Linewidth',2);
% plot(r3,'o-','Color',color_save_new(3),'MarkerEdgeColor',color_save_new(3), 'MarkerFaceColor',color_save_new(3),'Linewidth',2);
plot(r, 'o-','Color',color_save_new(5),'MarkerEdgeColor',color_save_new(5), 'MarkerFaceColor',color_save_new(5),'Linewidth',2);
xlabel('m','Fontsize',16);
ylabel('$\triangle_m (\%)$','Interpreter','Latex','Fontsize',16);
set(gca,'Xtick',[1,2,3,4,5,6,7],'Xticklabel',{'I','II','III','IV', 'V', 'VI', 'VII'});

% Legend Plotting
% ylim([0 10]);
% set(gca,'Fontsize',16);
% set(gca,'Fontname','Times');
% print('-depsc','error_check.eps');
% hh = legend('Train:Osc I','Train:Osc II','Aggregate','Location','Southoutside');
% legend boxoff;
% set(hh,'Fontname','Times');
% print('-depsc','error_check1.eps');

plotting_adjacency_new1(abs(Areg),7,1);
plotting_adjacency_new1(angle(Areg),7,2);


function [deg] = degree(A,type)
    % type: "in" or "out"

    if type == 'in'
        type = 2;
    elseif type == 'out'
        type = 1;
    else
        fprintf('Please specify either "in" or "out" as type! \n')
    end

    deg = sum(sum(abs(A),type));

end
