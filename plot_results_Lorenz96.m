clc,clear;

m1 = matfile('TDRC_Lorenz96_n5_N1200_delay2');
m2 = matfile('TDRC_Lorenz96_n5_N800_delay3');
m3 = matfile('TDRC_Lorenz96_n5_N700_delay5');

figure('position', [500,400,600,500])
p1 = plot(m1.ly_history(end,:),'b+','linewidth',2);hold on;
set(get(p1,'parent'),'linewidth',1.9)
p2 = plot(m2.ly_history(end,:),'rx','linewidth',2);hold on
p3 = plot(m3.ly_history(end,:),'mo','linewidth',2);hold on
lgd = legend([p1,p2,p3], {'$N=1200,\tau=2\Delta t$', '$N=800,\tau=3\Delta t$', '$N=700,\tau=5\Delta t$'}, 'interpreter','latex','fontsize', 20, 'box', 'off');  
lgd.ItemTokenSize = [60,20];  % set the legend length
p4 = plot([0, 21], [0,0], 'k--', 'linewidth', 1.8);hold off
set(get(get(p4, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
set(gca,'Position',[0.17,0.17,0.75,0.75],'fontsize',22)
xlabel('$j$', 'interpreter', 'latex', 'fontsize', 24);
ylabel('$\Lambda_j$', 'interpreter', 'latex', 'fontsize', 24)
xlim([0, 21])
ylim([-17, 17])


figure('name', 'prediction', 'position', [400,250,900,700])
subplot(511)
p1=plot(m3.t_pred, m3.u_target(1,:), 'r-', 'linewidth', 2);hold on
plot(m3.t_pred, m3.prediction(1,:), 'b-','linewidth', 2)
lgd = legend({'True$\quad$', 'prediction'}, 'interpreter','latex','fontsize', 22, 'box', 'off');  
lgd.ItemTokenSize = [60,20];  % set the legend length
set(lgd, 'Orientation', 'horizon')
set(get(p1,'parent'),'linewidth',1.9)  % set linewidth of axes
set(gca,'Position',[0.1,0.76,0.88,0.16],'fontsize',22)
set(gca,'XTicklabel',[]);
ylabel('$x_1$','interpreter','latex','fontsize',22)
ylim([-8, 12])
subplot(512)
p1=plot(m3.t_pred, m3.u_target(2,:), 'r-', 'linewidth', 2);hold on
plot(m3.t_pred, m3.prediction(2,:), 'b-', 'linewidth', 2)
set(get(p1,'parent'),'linewidth',1.9)
set(gca,'Position',[0.1,0.6,0.88,0.16],'fontsize',22)
set(gca,'XTicklabel',[]);
ylabel('$x_2$','interpreter','latex','fontsize',22)
ylim([-8, 12])
subplot(513)
p1=plot(m3.t_pred, m3.u_target(3,:), 'r-', 'linewidth', 2);hold on
plot(m3.t_pred, m3.prediction(3,:), 'b-', 'linewidth', 2)
set(get(p1,'parent'),'linewidth',1.9)
set(gca,'Position',[0.1,0.44,0.88,0.16],'fontsize',22)
set(gca,'XTicklabel',[]);
ylabel('$x_3$','interpreter','latex','fontsize',22)
ylim([-8, 12])
subplot(514)
p1=plot(m3.t_pred, m3.u_target(4,:), 'r-', 'linewidth', 2);hold on
plot(m3.t_pred, m3.prediction(4,:), 'b-', 'linewidth', 2)
set(get(p1,'parent'),'linewidth',1.9)
set(gca,'Position',[0.1,0.28,0.88,0.16],'fontsize',22)
set(gca,'XTicklabel',[]);
ylabel('$x_4$','interpreter','latex','fontsize',22)
subplot(515)
p1=plot(m3.t_pred, m3.u_target(5,:), 'r-', 'linewidth', 2);hold on
plot(m3.t_pred, m3.prediction(5,:), 'b-', 'linewidth', 2)
set(get(p1,'parent'),'linewidth',1.9)
set(gca,'Position',[0.1,0.12,0.88,0.16],'fontsize',22)
ylabel('$x_5$','interpreter','latex','fontsize',22)
xlabel('$\Lambda_1 t$','interpreter','latex','fontsize',22)
ylim([-8, 12])



[~, input_dim] = size(m1.u_target);
error1 = sqrt(sum((m1.u_target - m1.prediction).^2, 1)) / input_dim;
error2 = sqrt(sum((m2.u_target - m2.prediction).^2, 1)) / input_dim;
error3 = sqrt(sum((m3.u_target - m3.prediction).^2, 1)) / input_dim;

figure('position', [800,300,800,500])
p1 = plot(m1.t_pred, error1, 'b-','linewidth',2);hold on;
set(get(p1,'parent'),'linewidth',1.9)
p2 = plot(m2.t_pred, error2 ,'r-','linewidth',2);hold on
p3 = plot(m3.t_pred, error3 ,'m-','linewidth',2);hold on
lgd = legend([p1,p2,p3], {'$N=1200,\tau=2\Delta t$', ...
                          '$N=800,\tau=3\Delta t$', ...
                          '$N=700,\tau=5\Delta t$'}, ...
                          'interpreter','latex','fontsize', 20, ...
                          'location', 'northwest', 'box', 'off');  
lgd.ItemTokenSize = [40,10];  % set the legend length
set(gca,'Position',[0.19,0.17,0.75,0.75],'fontsize',22)
xlabel('$\Lambda_1 t$', 'interpreter', 'latex', 'fontsize', 24);
ylabel('NMSE', 'interpreter', 'latex', 'fontsize', 24)
xlim([0, 10])
ylim([-0.001, 0.014])












