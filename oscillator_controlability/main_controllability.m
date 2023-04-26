clear;clc;close all;
 
%% load Laplacian matrix
% ------------------------------------------------------------------
load aggrigate_L.mat %Lr_final.mat
A  = Lreg;

%% Initial parameters
% ------------------------------------------------------------------
np = 101; % number or runs
rho = logspace(-2,1,np); % "sigma" in paper

ocls = 7;

% LQR stuff 
Q = eye(ocls,ocls);
R = 1;
N = zeros(ocls,1);

nosc = ocls;

% color map stuff
c1 = hot(115);
c1 = c1(1:101,:);
c1 = flipud(c1);

%% Allocate Space
% ------------------------------------------------------------------
lam1    = zeros(nosc,np);
lam2    = zeros(nosc,np);
lam3    = zeros(nosc,np);
% lam4    = zeros(nosc,np);
% lam5    = zeros(nosc,np);
% lam6    = zeros(nosc,np);
% lam7    = zeros(nosc,np);
lam12   = zeros(nosc,np);
lam13   = zeros(nosc,np);
lam23   = zeros(nosc,np);
lam123 = zeros(nosc,np);
lam14    = zeros(nosc,np);
lam24    = zeros(nosc,np);
lam15    = zeros(nosc,np);
lam25    = zeros(nosc,np);

% --- structure modes ---

B  = [1 0 0 0 0 0 0]'; % input vector
for i = 1:length(rho)
    K = lqr(A,B,Q,R*rho(i),N);
    lam1(:,i) = eig(A-B*K);
end
B  = [0 1 0 0 0 0 0]';
for i = 1:length(rho)
    K = lqr(A,B,Q,R*rho(i),N);
    lam2(:,i) = eig(A-B*K);
end
B  = [0 0 1 0 0 0 0]';
for i = 1:length(rho)
    K = lqr(A,B,Q,R*rho(i),N);
    lam3(:,i) = eig(A-B*K);
end

% --- fluid modes ---

% B  = [0 0 0 1 0 0 0]';
% for i = 1:length(rho)
%     K = lqr(A,B,Q,R*rho(i),N);
%     lam4(:,i) = eig(A-B*K);
% end
% B  = [0 0 0 0 1 0 0]';
% for i = 1:length(rho)
%     K = lqr(A,B,Q,R*rho(i),N);
%     lam5(:,i) = eig(A-B*K);
% end
% B  = [0 0 0 0 0 1 0]';
% for i = 1:length(rho)
%     K = lqr(A,B,Q,R*rho(i),N);
%     lam6(:,i) = eig(A-B*K);
% end
% B  = [0 0 0 0 0 0 1]';
% for i = 1:length(rho)
%     K = lqr(A,B,Q,R*rho(i),N);
%     lam7(:,i) = eig(A-B*K);
% end

% --- multiple strucutre modes ---

B  = [1 1 0 0 0 0 0]';
for i = 1:length(rho)
    K = lqr(A,B,Q,R*rho(i),N);
    lam12(:,i) = eig(A-B*K);
end
B  = [1 0 1 0 0 0 0]';
for i = 1:length(rho)
    K = lqr(A,B,Q,R*rho(i),N);
    lam13(:,i) = eig(A-B*K);
end
B  = [0 1 1 0 0 0 0]';
for i = 1:length(rho)
    K = lqr(A,B,Q,R*rho(i),N);
    lam23(:,i) = eig(A-B*K);
end
B  = [1 1 1 0 0 0 0]';
for i = 1:length(rho)
    K = lqr(A,B,Q,R*rho(i),N);
    lam123(:,i) = eig(A-B*K);
end

% --- Fluid and Structure ---

B  = [1 0 0 1 0 0 0]';
for i = 1:length(rho)
    K = lqr(A,B,Q,R*rho(i),N);
    lam14(:,i) = eig(A-B*K);
end
B  = [1 0 0 0 1 0 0]';
for i = 1:length(rho)
    K = lqr(A,B,Q,R*rho(i),N);
    lam15(:,i) = eig(A-B*K);
end
B  = [0 1 0 1 0 0 0]';
for i = 1:length(rho)
    K = lqr(A,B,Q,R*rho(i),N);
    lam24(:,i) = eig(A-B*K);
end
B  = [0 1 0 0 1 0 0]';
for i = 1:length(rho)
    K = lqr(A,B,Q,R*rho(i),N);
    lam25(:,i) = eig(A-B*K);
end

fs = 18;
lw = 2;

%% Plotting
% ------------------------------------------------------------------

% --- structure ---

figure;
subplot(221);
plot(eig(A),'*','Color','k','Markersize',12);hold on;
for i = 1:length(rho)
    plot(lam1(:,i),'o','Color',c1(i,:));hold on;
end
xlim([-12 0]);ylim([-10 10]);
caxis(log([rho(1) rho(length(rho))]));
%colorbar('FontSize',11,'YTick',[log(0.1) log(1) log(10) log(100) log(1000)],'YTickLabel',[0.1 1 10 100 1000]);
%title('B = [1 0 0 0]');
xlabel('$\Re(\lambda)$','Interpreter','Latex','Fontsize',fs);
ylabel('$\Im(\lambda)$','Interpreter','Latex','Fontsize',fs);
plotTickLatex2D('Fontsize',fs);
box off;grid on;
print('-depsc','eig_move1.eps');

figure;
subplot(221);
plot(eig(A),'*','Color','k','Markersize',12);hold on;
for i = 1:length(rho)
    plot(lam2(:,i),'o','Color',c1(i,:));hold on;
end
xlim([-12 0]);ylim([-10 10]);
caxis(log([rho(1) rho(length(rho))]));
%colorbar('FontSize',11,'YTick',[log(0.1) log(1) log(10) log(100) log(1000)],'YTickLabel',[0.1 1 10 100 1000]);
%title('B = [0 1 0 0]');
xlabel('$\Re(\lambda)$','Interpreter','Latex','Fontsize',fs);
ylabel('$\Im(\lambda)$','Interpreter','Latex','Fontsize',fs);
plotTickLatex2D('Fontsize',fs);
box off;grid on;
print('-depsc','eig_move2.eps');

figure;
subplot(221);
plot(eig(A),'*','Color','k','Markersize',12);hold on;
for i = 1:length(rho)
    plot(lam3(:,i),'o','Color',c1(i,:));hold on;
end
xlim([-12 0]);ylim([-10 10]);
caxis(log([rho(1) rho(length(rho))]));
%colorbar('FontSize',11,'YTick',[log(0.1) log(1) log(10) log(100) log(1000)],'YTickLabel',[0.1 1 10 100 1000]);
%title('B = [0 0 1 0]');
xlabel('$\Re(\lambda)$','Interpreter','Latex','Fontsize',fs);
ylabel('$\Im(\lambda)$','Interpreter','Latex','Fontsize',fs);
plotTickLatex2D('Fontsize',fs);
box off;grid on;
print('-depsc','eig_move3.eps');

% --- Fluid modes ---

% figure;
% subplot(221);
% plot(eig(A),'*','Color','k','Markersize',12);hold on;
% for i = 1:length(rho)
%     plot(lam4(:,i),'o','Color',c1(i,:));hold on;
% end
% xlim([-12 0]);ylim([-10 10]);
% %caxis(log([rho(1) rho(length(rho))]));
% %colorbar('FontSize',11,'XTick',[1 26 51 76 101],'YTickLabel',[0.1 1 10 100 1000]);
% %title('B = [0 0 0 1]');
% xlabel('$\Re(\lambda)$','Interpreter','Latex','Fontsize',fs);
% ylabel('$\Im(\lambda)$','Interpreter','Latex','Fontsize',fs);
% plotTickLatex2D('Fontsize',fs);
% box off;grid on;
% print('-depsc','eig_move4.eps');
% 
% figure;
% subplot(221);
% plot(eig(A),'*','Color','k','Markersize',12);hold on;
% for i = 1:length(rho)
%     plot(lam5(:,i),'o','Color',c1(i,:));hold on;
% end
% xlim([-12 0]);ylim([-10 10]);
% %caxis(log([rho(1) rho(length(rho))]));
% %colorbar('FontSize',11,'XTick',[1 26 51 76 101],'YTickLabel',[0.1 1 10 100 1000]);
% %title('B = [0 0 0 1]');
% xlabel('$\Re(\lambda)$','Interpreter','Latex','Fontsize',fs);
% ylabel('$\Im(\lambda)$','Interpreter','Latex','Fontsize',fs);
% plotTickLatex2D('Fontsize',fs);
% box off;grid on;
% print('-depsc','eig_move4.eps');
% 
% figure;
% subplot(221);
% plot(eig(A),'*','Color','k','Markersize',12);hold on;
% for i = 1:length(rho)
%     plot(lam6(:,i),'o','Color',c1(i,:));hold on;
% end
% xlim([-12 0]);ylim([-10 10]);
% %caxis(log([rho(1) rho(length(rho))]));
% %colorbar('FontSize',11,'XTick',[1 26 51 76 101],'YTickLabel',[0.1 1 10 100 1000]);
% %title('B = [0 0 0 1]');
% xlabel('$\Re(\lambda)$','Interpreter','Latex','Fontsize',fs);
% ylabel('$\Im(\lambda)$','Interpreter','Latex','Fontsize',fs);
% plotTickLatex2D('Fontsize',fs);
% box off;grid on;
% print('-depsc','eig_move4.eps');
% 
% figure;
% subplot(221);
% plot(eig(A),'*','Color','k','Markersize',12);hold on;
% for i = 1:length(rho)
%     plot(lam7(:,i),'o','Color',c1(i,:));hold on;
% end
% xlim([-12 0]);ylim([-10 10]);
% %caxis(log([rho(1) rho(length(rho))]));
% %colorbar('FontSize',11,'XTick',[1 26 51 76 101],'YTickLabel',[0.1 1 10 100 1000]);
% %title('B = [0 0 0 1]');
% xlabel('$\Re(\lambda)$','Interpreter','Latex','Fontsize',fs);
% ylabel('$\Im(\lambda)$','Interpreter','Latex','Fontsize',fs);
% plotTickLatex2D('Fontsize',fs);
% box off;grid on;
% print('-depsc','eig_move4.eps');

% --- multiple structure modes ---

figure;
subplot(221);
plot(eig(A),'*','Color','k','Markersize',12);hold on;
for i = 1:length(rho)
    plot(lam12(:,i),'o','Color',c1(i,:));hold on;
end
xlim([-12 0]);ylim([-10 10]);
caxis(log([rho(1) rho(length(rho))]));
xlabel('$\Re(\lambda)$','Interpreter','Latex','Fontsize',fs);
ylabel('$\Im(\lambda)$','Interpreter','Latex','Fontsize',fs);
plotTickLatex2D('Fontsize',fs);
box off;grid on;
print('-depsc','eig_move12.eps');

figure;
subplot(221);
plot(eig(A),'*','Color','k','Markersize',12);hold on;
for i = 1:length(rho)
    plot(lam13(:,i),'o','Color',c1(i,:));hold on;
end
xlim([-12 0]);ylim([-10 10]);
caxis(log([rho(1) rho(length(rho))]));
xlabel('$\Re(\lambda)$','Interpreter','Latex','Fontsize',fs);
ylabel('$\Im(\lambda)$','Interpreter','Latex','Fontsize',fs);
plotTickLatex2D('Fontsize',fs);
box off;grid on;
print('-depsc','eig_move13.eps');

figure;
subplot(221);
plot(eig(A),'*','Color','k','Markersize',12);hold on;
for i = 1:length(rho)
    plot(lam23(:,i),'o','Color',c1(i,:));hold on;
end
xlim([-12 0]);ylim([-10 10]);
caxis(log([rho(1) rho(length(rho))]));
xlabel('$\Re(\lambda)$','Interpreter','Latex','Fontsize',fs);
ylabel('$\Im(\lambda)$','Interpreter','Latex','Fontsize',fs);
plotTickLatex2D('Fontsize',fs);
box off;grid on;
print('-depsc','eig_move14.eps');

figure;
subplot(221);
plot(eig(A),'*','Color','k','Markersize',12);hold on;
for i = 1:length(rho)
    plot(lam123(:,i),'o','Color',c1(i,:));hold on;
end
xlim([-12 0]);ylim([-10 10]);
caxis(log([rho(1) rho(length(rho))]));
xlabel('$\Re(\lambda)$','Interpreter','Latex','Fontsize',fs);
ylabel('$\Im(\lambda)$','Interpreter','Latex','Fontsize',fs);
plotTickLatex2D('Fontsize',fs);
box off;grid on;
print('-depsc','eig_move1234.eps');

% --- Fluid and Structure ---

figure;
subplot(221);
plot(eig(A),'*','Color','k','Markersize',12);hold on;
for i = 1:length(rho)
    plot(lam14(:,i),'o','Color',c1(i,:));hold on;
end
xlim([-12 0]);ylim([-10 10]);
caxis(log([rho(1) rho(length(rho))]));
xlabel('$\Re(\lambda)$','Interpreter','Latex','Fontsize',fs);
ylabel('$\Im(\lambda)$','Interpreter','Latex','Fontsize',fs);
plotTickLatex2D('Fontsize',fs);
box off;grid on;
print('-depsc','eig_move1234.eps');

figure;
subplot(221);
plot(eig(A),'*','Color','k','Markersize',12);hold on;
for i = 1:length(rho)
    plot(lam15(:,i),'o','Color',c1(i,:));hold on;
end
xlim([-12 0]);ylim([-10 10]);
caxis(log([rho(1) rho(length(rho))]));
xlabel('$\Re(\lambda)$','Interpreter','Latex','Fontsize',fs);
ylabel('$\Im(\lambda)$','Interpreter','Latex','Fontsize',fs);
plotTickLatex2D('Fontsize',fs);
box off;grid on;
print('-depsc','eig_move1234.eps');

figure;
subplot(221);
plot(eig(A),'*','Color','k','Markersize',12);hold on;
for i = 1:length(rho)
    plot(lam24(:,i),'o','Color',c1(i,:));hold on;
end
xlim([-12 0]);ylim([-10 10]);
caxis(log([rho(1) rho(length(rho))]));
xlabel('$\Re(\lambda)$','Interpreter','Latex','Fontsize',fs);
ylabel('$\Im(\lambda)$','Interpreter','Latex','Fontsize',fs);
plotTickLatex2D('Fontsize',fs);
box off;grid on;
print('-depsc','eig_move1234.eps');

figure;
subplot(221);
plot(eig(A),'*','Color','k','Markersize',12);hold on;
for i = 1:length(rho)
    plot(lam25(:,i),'o','Color',c1(i,:));hold on;
end
xlim([-12 0]);ylim([-10 10]);
caxis(log([rho(1) rho(length(rho))]));
xlabel('$\Re(\lambda)$','Interpreter','Latex','Fontsize',fs);
ylabel('$\Im(\lambda)$','Interpreter','Latex','Fontsize',fs);
plotTickLatex2D('Fontsize',fs);
box off;grid on;
print('-depsc','eig_move1234.eps');


%% Maximum Pole Movement
% ------------------------------------------------------------------

% lam1_max     = zeros(length(rho),1);
% lam12_max    = zeros(length(rho),1);
% lam13_max    = zeros(length(rho),1);
% lam23_max    = zeros(length(rho),1);
% lam1234_max  = zeros(length(rho),1);
% lami1_max    = zeros(length(rho),1);
% lami12_max   = zeros(length(rho),1);
% lami13_max   = zeros(length(rho),1);
% lami23_max   = zeros(length(rho),1);
% lami1234_max = zeros(length(rho),1);
%
% idx1    = zeros(length(rho),1);
% idx12   = zeros(length(rho),1);
% idx13   = zeros(length(rho),1);
% idx23   = zeros(length(rho),1);
% idx1234 = zeros(length(rho),1);
% 
% for i = 1:length(rho)
%     [lam1_max(i), idx1(i)] = max(real(lam1(:,i)));
%     lami1_max(i) = imag(lam1(idx1(i),i));
% 
%     [lam12_max(i), idx12(i)] = max(real(lam12(:,i)));
%     lami12_max(i) = imag(lam12(idx12(i),i));
% 
%     [lam13_max(i), idx13(i)] = max(real(lam13(:,i)));
%     lami13_max(i) = imag(lam13(idx13(i),i));
% 
%     [lam23_max(i), idx23(i)] = max(real(lam23(:,i)));
%     lami23_max(i) = imag(lam23(idx23(i),i));
% 
%     [lam1234_max(i), idx1234(i)] = max(real(lam1234(:,i)));
%     lami1234_max(i) = imag(lam1234(idx1234(i),i));
% end

%%  Maximum Pole Movement Plotting
% ------------------------------------------------------------------
end_mode = 4;

figure;subplot(211); 
plot(rho, max(real(lam12(1:end_mode,:))),'r-','Linewidth',lw); hold on;
plot(rho, max(real(lam13(1:end_mode,:))),'k-','Linewidth',lw); hold on;
plot(rho, max(real(lam23(1:end_mode,:))),'b-','Linewidth',lw); hold on;
plot(rho, max(real(lam123(1:end_mode,:))),'--','color', [0.1843 0.7412 0.3529],'Linewidth',lw); hold on;
set(gca,'Xscale','log');
xlim([10^-1 10^2]);
set(gca,'Xtick',[0.1 1 10 100 1000]);
%hc = legend('1','12','13','14',2);set(hc,'Fontsize',14);
xlabel('$\sigma$','Interpreter','Latex','Fontsize',14);
ylabel('$max(\Re(\lambda))$','Interpreter','Latex','Fontsize',14);
legend boxoff;
legend('[1 1 0 0 0 0 0]', '[1 0 1 0 0 0 0]', '[0 1 1 0 0 0 0]', '[1 1 1 1 1 1 1]')
plotTickLatex2D('Fontsize',14);
box off;grid on;set(gca,'XMinorGrid','Off');
print('-depsc','eig_max.eps');


% % % 
% % % figure;subplot(311); 
% % % %plot(max(real(eig(A))),'*','Color','k');hold on;
% % % plot(rho,lam12_max,'r-','Linewidth',lw);hold on;
% % % plot(rho,lam13_max,'k-','Linewidth',lw);hold on;
% % % plot(rho,lam23_max,'b-','Linewidth',lw);hold on;
% % % plot(rho,lam1234_max,'g-','Linewidth',lw);hold on;
% % % %plot(rho,lam1234_max,'m-','Linewidth',lw);hold on;
% % % %plot(rho,lam1234_max,'');hold on;
% % % set(gca,'Xscale','log');
% % % xlim([10^-1 10^2]);
% % % set(gca,'Xtick',[0.1 1 10 100 1000]);
% % % %hc = legend('1','12','13','14',2);set(hc,'Fontsize',14);
% % % xlabel('$\sigma$','Interpreter','Latex','Fontsize',14);
% % % ylabel('$max(\Re(\lambda))$','Interpreter','Latex','Fontsize',14);
% % % legend boxoff;
% % % plotTickLatex2D('Fontsize',14);
% % % box off;grid on;set(gca,'XMinorGrid','Off');
% % % print('-depsc','eig_max.eps');

% figure;subplot(311); 
% %plot(max(real(eig(A))),'*','Color','k');hold on;
% plot(rho,lam12_max,'r-','Linewidth',lw);hold on;
% plot(rho,lam13_max,'k-','Linewidth',lw);hold on;
% plot(rho,lam23_max,'b-','Linewidth',lw);hold on;
% plot(rho,lam1234_max,'Color',[0.4660, 0.6740, 0.1880],'Linewidth',lw);hold on;
% %plot(rho,lam1234_max,'m-','Linewidth',lw);hold on;
% %plot(rho,lam1234_max,'');hold on;
% set(gca,'Xscale','log');
% xlim([10^-1 10^3]);
% set(gca,'Xtick',[0.1 1 10 100 1000]);
% % % % % % % hc = legend('12','13','14','1234',2);set(hc,'Fontsize',14);
% xlabel('\rho','Fontsize',14);ylabel('max(\lambda_r)','Fontsize',14);
% %legend boxoff;
% plotTickLatex2D('Fontsize',14);
% box off;grid on;set(gca,'XMinorGrid','Off');
% print('-depsc','eig_max_leg.eps');
% 
% 
% figure;subplot(311); 
% plot(rho,lami1_max,'r-','Linewidth',lw);hold on;
% plot(rho,lami12_max,'k-','Linewidth',lw);hold on;
% plot(rho,lami13_max,'b-','Linewidth',lw);hold on;
% plot(rho,lami23_max,'Color',[0.4660, 0.6740, 0.1880],'Linewidth',lw);hold on;
% %plot(rho,lam1234_max,'Color',[0, 0.6740, 0],'Linewidth',lw);hold on;
% %plot(rho,lam1234_max,'');hold on;
% set(gca,'Xscale','log');
% xlim([10^-1 10^3]);
% set(gca,'Xtick',[0.1 1 10 100 1000]);
% %hc = legend('1','12','13','14',2);set(hc,'Fontsize',14);
% xlabel('\rho','Fontsize',14);ylabel('max(\lambda_r)','Fontsize',14);
% legend boxoff;
% plotTickLatex2D('Fontsize',14);
% box off;grid on;set(gca,'XMinorGrid','Off');
% print('-depsc','eig_max_imag.eps');

% ------------------------------------------------------------------
% ------------------------------------------------------------------
% ------------------------------------------------------------------
% ------------------------------------------------------------------
% ------------------------------------------------------------------


% % % % lam1_l2 = zeros(length(rho),1);
% % % % lam2_l2 = zeros(length(rho),1);
% % % % lam3_l2 = zeros(length(rho),1);
% % % % lam4_l2 = zeros(length(rho),1);
% % % % idx1_1 = 2*ones(29,1);
% % % % idx1_2 = 1*ones(57,1);
% % % % idx1_3 = 2*ones(15,1);
% % % % idx1_l2 = [idx1_1;idx1_2;idx1_3];
% % % % for i = 1:length(rho)
% % % %     lam1_l2(i) = real(lam1(idx1_l2(i),i));
% % % % end
% % % % idx2_1 = 1*ones(17,1);
% % % % idx2_2 = 2*ones(84,1);
% % % % idx2_l2 = [idx2_1;idx2_2];
% % % % for i = 1:length(rho)
% % % %     lam2_l2(i) = real(lam2(idx2_l2(i),i));
% % % % end
% % % % idx3_1 = 3*ones(22,1);
% % % % idx3_2 = 2*ones(79,1);
% % % % idx3_l2 = [idx3_1;idx3_2];
% % % % for i = 1:length(rho)
% % % %     lam3_l2(i) = real(lam3(idx3_l2(i),i));
% % % % end
% % % % idx4_1 = 3*ones(48,1);
% % % % idx4_2 = 2*ones(53,1);
% % % % idx4_l2 = [idx4_1;idx4_2];
% % % % for i = 1:length(rho)
% % % %     lam4_l2(i) = real(lam4(idx4_l2(i),i));
% % % % end
% % % % 
% % % % 
% % % % figure;subplot(311); 
% % % % %plot(max(real(eig(A))),'*','Color','k');hold on;
% % % % plot(rho,lam1_l2,'r-','Linewidth',lw);hold on;
% % % % plot(rho,lam2_l2,'k-','Linewidth',lw);hold on;
% % % % plot(rho,lam3_l2,'b-','Linewidth',lw);hold on;
% % % % plot(rho,lam4_l2,'Color',[0.4660, 0.6740, 0.1880],'Linewidth',lw);hold on;
% % % % %plot(rho,lam1234_max,'');hold on;
% % % % set(gca,'Xscale','log');
% % % % xlim([10^-1 10^3]);
% % % % set(gca,'Xtick',[0.1 1 10 100 1000]);
% % % % %hc = legend('1','12','13','14','1234',2);set(hc,'Fontsize',14);
% % % % xlabel('\rho','Fontsize',14);ylabel('(\lambda_r)_{II}','Fontsize',14);
% % % % %legend boxoff;
% % % % plotTickLatex2D('Fontsize',14);
% % % % box off;grid on;
% % % % print('-depsc','eig_2.eps');
% % % % 
% % % % figure;subplot(311); 
% % % % %plot(max(real(eig(A))),'*','Color','k');hold on;
% % % % plot(rho,lam1_l2,'r-','Linewidth',lw);hold on;
% % % % plot(rho,lam2_l2,'k-','Linewidth',lw);hold on;
% % % % plot(rho,lam3_l2,'b-','Linewidth',lw);hold on;
% % % % plot(rho,lam4_l2,'Color',[0.4660, 0.6740, 0.1880],'Linewidth',lw);hold on;
% % % % %plot(rho,lam1234_max,'');hold on;
% % % % set(gca,'Xscale','log');
% % % % xlim([10^-1 10^3]);
% % % % set(gca,'Xtick',[0.1 1 10 100 1000]);
% % % % % % % % % % hc = legend('1','2','3','4',2);set(hc,'Fontsize',14);
% % % % xlabel('\rho','Fontsize',14);ylabel('\lambda_{II}','Fontsize',14);
% % % % %legend boxoff;
% % % % plotTickLatex2D('Fontsize',14);
% % % % box off;grid on;
% % % % print('-depsc','eig_2_leg.eps');

