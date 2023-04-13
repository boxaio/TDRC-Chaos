clc,clear;


% results of Lorenz96 (n = 5)

dt = 0.02;
set_average_degree = 3;
approx_reservoir_size = 700;
rho = 1.2;
num_delay = 5;
gama = 0.1;
beta = 1e-6;
c = 0.1;
num_inter_step = 20;

compute_lyap = 0;
num_lyap = 20;


[t_pred, u_target, prediction] ...
           = TDRC_Lorenz96_n5(dt, set_average_degree, approx_reservoir_size,...
                              rho, num_delay, gama, beta, c, num_inter_step, ...
                              compute_lyap, num_lyap);
                          
figure('name', 'prediction', 'position', [200,50,1100,1000])
subplot(511)
p1=plot(t_pred(1:1000), u_target(1,1:1000), 'r-', 'linewidth', 2);hold on
plot(t_pred(1:1000), prediction(1,1:1000), 'b-','linewidth', 2)
lgd = legend({'True$\quad$', 'prediction'}, 'interpreter','latex','fontsize', 22, 'box', 'off');  
lgd.ItemTokenSize = [60,20];  % set the legend length
set(lgd, 'Orientation', 'horizon')
set(get(p1,'parent'),'linewidth',1.9)  % set linewidth of axes
set(gca,'Position',[0.08,0.78,0.9,0.17],'fontsize',22)
set(gca,'XTicklabel',[]);
ylabel('$x_1$','interpreter','latex','fontsize',22)
ylim([-8, 12])
subplot(512)
p1=plot(t_pred(1:1000), u_target(2,1:1000), 'r-', 'linewidth', 2);hold on
plot(t_pred(1:1000), prediction(2,1:1000), 'b-', 'linewidth', 2)
set(get(p1,'parent'),'linewidth',1.9)
set(gca,'Position',[0.08,0.61,0.9,0.17],'fontsize',22)
set(gca,'XTicklabel',[]);
ylabel('$x_2$','interpreter','latex','fontsize',22)
ylim([-8, 12])
subplot(513)
p1=plot(t_pred(1:1000), u_target(3,1:1000), 'r-', 'linewidth', 2);hold on
plot(t_pred(1:1000), prediction(3,1:1000), 'b-', 'linewidth', 2)
set(get(p1,'parent'),'linewidth',1.9)
set(gca,'Position',[0.08,0.44,0.9,0.17],'fontsize',22)
set(gca,'XTicklabel',[]);
ylabel('$x_3$','interpreter','latex','fontsize',22)
ylim([-8, 12])
subplot(514)
p1=plot(t_pred(1:1000), u_target(4,1:1000), 'r-', 'linewidth', 2);hold on
plot(t_pred(1:1000), prediction(4,1:1000), 'b-', 'linewidth', 2)
set(get(p1,'parent'),'linewidth',1.9)
set(gca,'Position',[0.08,0.27,0.9,0.17],'fontsize',22)
set(gca,'XTicklabel',[]);
ylabel('$x_4$','interpreter','latex','fontsize',22)
subplot(515)
p1=plot(t_pred(1:1000), u_target(5,1:1000), 'r-', 'linewidth', 2);hold on
plot(t_pred(1:1000), prediction(5,1:1000), 'b-', 'linewidth', 2)
set(get(p1,'parent'),'linewidth',1.9)
set(gca,'Position',[0.08,0.1,0.9,0.17],'fontsize',22)
ylabel('$x_5$','interpreter','latex','fontsize',22)
xlabel('$\Lambda_1 t$','interpreter','latex','fontsize',22)
ylim([-8, 12])


compute_lyap = 1;

[t_pred, u_target, prediction, ly_next, ly_history] ...
            = TDRC_Lorenz96_n5(dt, set_average_degree, approx_reservoir_size,...
                               rho, num_delay, gama, beta, c, num_inter_step, ...
                               compute_lyap, num_lyap);

figure('name', 'LZ96_n5 Lyapunov', 'position', [500,400,800,700])
plot(ly_history(:,1), 'linewidth', 2);hold on
plot(ly_history(:,2), 'linewidth', 2);hold on
plot(ly_history(:,3), 'linewidth', 2);hold on
plot(ly_history(:,4), 'linewidth', 2);hold on
plot(ly_history(:,5), 'linewidth', 2);hold on
set(gca,'Position',[0.1,0.11,0.85,0.85],'fontsize',20)
xlabel('Iterations', 'fontsize', 22);
ylabel('$\Lambda_j$', 'interpreter', 'latex', 'fontsize', 22)

save(['TDRC_Lorenz96_n5', '_N', num2str(approx_reservoir_size), '_delay', num2str(num_delay), '.mat'],...
      'dt','approx_reservoir_size','num_delay','rho','gama','beta', 'c','num_inter_step', 'num_lyap', 'u_target','t_pred', 'prediction', 'ly_history')
  
% figure('position', [200,400,800,700])
% p1 = plot(ly_next, 'b+', 'linewidth', 2);hold on
% set(get(p1, 'parent'), 'linewidth', 1.9)
% plot([0,num_lyap], [0,0], 'k--', 'linewidth', 1.8);
% set(gca,'Position',[0.11,0.11,0.85,0.85],'fontsize',20)
% xlabel('$j$', 'interpreter', 'latex', 'fontsize', 22);
% ylabel('$\Lambda_j$', 'interpreter', 'latex', 'fontsize', 22)
