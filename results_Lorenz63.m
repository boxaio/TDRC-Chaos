clc,clear;

% 4.386917861452845   3.008681540387677   0.381561841627546  -2.947404298063971  -3.092527438436751  -3.208329809722969
% results of Lorenz63

dt = 0.02;
set_average_degree = 3;
approx_reservoir_size = 100;
rho = 1.2;
num_delay = 20;
gama = 0.2;
beta = 1e-6;
c = 0.1;
num_inter_step = 20;

compute_lyap = 0;
num_lyap = 10;

[t_pred, u_target, prediction] ...
            = TDRC_Lorenz63(dt, set_average_degree, approx_reservoir_size,...
                            rho, num_delay, gama, beta, c, num_inter_step, ...
                            compute_lyap, num_lyap);

figure('name', 'prediction', 'position', [200,150,1100,700])
subplot(311)
p1=plot(t_pred, u_target(1,:), 'r-', 'linewidth', 2);hold on
plot(t_pred, prediction(1,:), 'b-','linewidth', 2)
lgd = legend({'True$\quad$', 'Prediction'}, 'interpreter','latex','fontsize', 22, 'box', 'off');  
lgd.ItemTokenSize = [60,20];  % set the legend length
set(lgd, 'Orientation', 'horizon')
set(get(p1,'parent'),'linewidth',1.9)  % set linewidth of axes
set(gca,'Position',[0.1,0.65,0.87,0.27],'fontsize',22)
set(gca,'XTicklabel',[]);
ylabel('$x_1$','interpreter','latex','fontsize',22)
xlim([0, 20])
ylim([-25, 25])
subplot(312)
p1=plot(t_pred, u_target(2,:), 'r-', 'linewidth', 2);hold on
plot(t_pred, prediction(2,:), 'b-', 'linewidth', 2)
set(get(p1,'parent'),'linewidth',1.9)
set(gca,'Position',[0.1,0.38,0.87,0.27],'fontsize',22)
set(gca,'XTicklabel',[]);
ylabel('$x_2$','interpreter','latex','fontsize',22)
xlim([0, 20])
ylim([-25, 25])
subplot(313)
p1=plot(t_pred, u_target(3,:), 'r-', 'linewidth', 2);hold on
plot(t_pred, prediction(3,:), 'b-', 'linewidth', 2)
set(get(p1,'parent'),'linewidth',1.9)
set(gca,'Position',[0.1,0.11,0.87,0.27],'fontsize',22)
xlabel('$\Lambda_1 t$','interpreter','latex','fontsize',22)
ylabel('$x_3$','interpreter','latex','fontsize',22)
xlim([0, 20])
ylim([0, 48])

compute_lyap = 1;

[t_pred, u_target, prediction, ly_next, ly_history] ...
            = TDRC_Lorenz63(dt, set_average_degree, approx_reservoir_size,...
                            rho, num_delay, gama, beta, c, num_inter_step, ...
                            compute_lyap, num_lyap);

figure('name', 'LZ63 Lyapunov', 'position', [500,400,800,700])
plot(ly_history(:,1), 'linewidth', 2);hold on
plot(ly_history(:,2), 'linewidth', 2);hold on
plot(ly_history(:,3), 'linewidth', 2);hold on
set(gca,'Position',[0.11,0.11,0.85,0.85],'fontsize',20)
xlabel('Iterations', 'fontsize', 22);
ylabel('$\Lambda_j$', 'interpreter', 'latex', 'fontsize', 22)

save(['TDRC_Lorenz63', '_N', num2str(approx_reservoir_size), '_delay', num2str(num_delay), '.mat'],...
      'dt','approx_reservoir_size','num_delay','rho','gama','beta', 'c','num_inter_step', 'num_lyap', 'u_target','t_pred', 'prediction', 'ly_history');
