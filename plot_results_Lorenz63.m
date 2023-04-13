clc,clear;

m1 = matfile('TDRC_Lorenz63_N200_delay5');
m2 = matfile('TDRC_Lorenz63_N150_delay10');
m3 = matfile('TDRC_Lorenz63_N100_delay20');

figure('position', [500,400,600,500])
p1 = plot(m1.ly_history(end,:),'b+','linewidth',2);hold on;
set(get(p1,'parent'),'linewidth',1.9)
p2 = plot(m2.ly_history(end,:),'rx','linewidth',2);hold on
p3 = plot(m3.ly_history(end,:),'mo','linewidth',2);hold on
lgd = legend([p1,p2,p3], {'$N=200,\tau=5\Delta t$', '$N=150,\tau=10\Delta t$', '$N=100,\tau=20\Delta t$'}, 'interpreter','latex','fontsize', 20, 'box', 'off');  
lgd.ItemTokenSize = [60,20];  % set the legend length
p4 = plot([0, 11], [0, 0], 'k--', 'linewidth', 1.8);hold off
set(get(get(p4, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
set(gca,'Position',[0.17,0.17,0.75,0.75],'fontsize',22)
xlabel('$j$', 'interpreter', 'latex', 'fontsize', 24);
ylabel('$\Lambda_j$', 'interpreter', 'latex', 'fontsize', 24)
xlim([0, 11])
ylim([-7, 7])


figure('name', 'prediction', 'position', [200,150,900,600])
subplot(311)
p1=plot(m3.t_pred, m3.u_target(1,:), 'r-', 'linewidth', 2);hold on
plot(m3.t_pred, m3.prediction(1,:), 'b-','linewidth', 2)
lgd = legend({'True$\quad$', 'Prediction'}, 'interpreter','latex','fontsize', 22, 'box', 'off');  
lgd.ItemTokenSize = [60,20];  % set the legend length
set(lgd, 'Orientation', 'horizon')
set(get(p1,'parent'),'linewidth',1.9)  % set linewidth of axes
set(gca,'Position',[0.1,0.67,0.87,0.27],'fontsize',22)
set(gca,'XTicklabel',[]);
ylabel('$x_1$','interpreter','latex','fontsize',22)
xlim([0, 20])
ylim([-25, 25])
subplot(312)
p1=plot(m3.t_pred, m3.u_target(2,:), 'r-', 'linewidth', 2);hold on
plot(m3.t_pred, m3.prediction(2,:), 'b-', 'linewidth', 2)
set(get(p1,'parent'),'linewidth',1.9)
set(gca,'Position',[0.1,0.4,0.87,0.27],'fontsize',22)
set(gca,'XTicklabel',[]);
ylabel('$x_2$','interpreter','latex','fontsize',22)
xlim([0, 20])
ylim([-25, 25])
subplot(313)
p1=plot(m3.t_pred, m3.u_target(3,:), 'r-', 'linewidth', 2);hold on
plot(m3.t_pred, m3.prediction(3,:), 'b-', 'linewidth', 2)
set(get(p1,'parent'),'linewidth',1.9)
set(gca,'Position',[0.1,0.13,0.87,0.27],'fontsize',22)
xlabel('$\Lambda_1 t$','interpreter','latex','fontsize',24)
ylabel('$x_3$','interpreter','latex','fontsize',24)
xlim([0, 20])
ylim([0, 48])

[~, input_dim] = size(m1.u_target);
error1 = sqrt(sum((m1.u_target - m1.prediction).^2, 1)) / input_dim;
error2 = sqrt(sum((m2.u_target - m2.prediction).^2, 1)) / input_dim;
error3 = sqrt(sum((m3.u_target - m3.prediction).^2, 1)) / input_dim;

figure('position', [800,300,800,500])
p1 = plot(m1.t_pred, error1, 'b-','linewidth',2);hold on;
set(get(p1,'parent'),'linewidth',1.9)
p2 = plot(m2.t_pred, error2 ,'r-','linewidth',2);hold on
p3 = plot(m3.t_pred, error3 ,'m-','linewidth',2);hold on
lgd = legend([p1,p2,p3], {'$N=200,\tau=5\Delta t$', ...
                          '$N=150,\tau=10\Delta t$', ...
                          '$N=100,\tau=20\Delta t$'}, ...
                          'interpreter','latex','fontsize', 20, ...
                          'location', 'northwest', 'box', 'off');  
lgd.ItemTokenSize = [40,10];  % set the legend length
set(gca,'Position',[0.19,0.17,0.75,0.75],'fontsize',22)
xlabel('$\Lambda_1 t$', 'interpreter', 'latex', 'fontsize', 24);
ylabel('NMSE', 'interpreter', 'latex', 'fontsize', 24)
xlim([0, 50])
ylim([-0.001, 0.018])














