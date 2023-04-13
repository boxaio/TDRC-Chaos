clc,clear;


% results of KS equation

m = matfile('KS_L22_mu007'); 
dt = m.dt;
L = m.L;
[~, input_dim] = size(m.ut);
set_average_degree = 6;
approx_reservoir_size = 800;
rho = 1.5;
num_delay = 5;
gama = 0.15;
beta = 1e-8;
c = 0.1;
num_inter_step = 20;

compute_lyap = 0;
num_lyap = 30;

[td, tt, xx, u_true, prediction] ...
            = TDRC_KS(m, set_average_degree, approx_reservoir_size,...
                      rho, num_delay, gama, beta, c, num_inter_step, ...
                      compute_lyap, num_lyap);
                  
error = u_true - prediction;

% Plot contour
figure('name', 'KS solution', 'position', [200,150,900,600])
subplot(311)
p1 = pcolor(tt, xx, u_true);
shading interp
colormap(slanCM('spectral', 128))
colorbar('Ticks', [-2,-1,0,1,2], 'position', [0.94, 0.2, 0.02, 0.6]);
set(get(p1, 'parent'), 'linewidth', 2)
title('True','fontsize',20)
set(gca,'XTicklabel',[]);
ylabel('$x$','interpreter','latex','fontsize',22)  
set(gca,'Position',[0.11,0.72,0.8,0.21],'fontsize',22)
xlim([-dt/10, td(end)+dt/10])
ylim([-L/300, L+L/300])
subplot(312)
p1 = pcolor(tt, xx, prediction);
shading interp
set(get(p1, 'parent'), 'linewidth', 2)
colormap(slanCM('spectral', 128))
title('Prediction','fontsize',20)
set(gca,'XTicklabel',[]);
ylabel('$x$','interpreter','latex','fontsize',22)  
set(gca,'Position',[0.11,0.43,0.8,0.21],'fontsize',22)
xlim([-dt/10, td(end)+dt/10])
ylim([-L/300, L+L/300])
subplot(313)
p1 = pcolor(tt, xx, error);
shading interp
set(get(p1, 'parent'), 'linewidth', 2)
colormap(slanCM('spectral', 128))
title('Error','fontsize',20)
ylabel('$x$','interpreter','latex','fontsize',22)  
set(gca,'Position',[0.11,0.13,0.8,0.21],'fontsize',22)
xlabel('$\Lambda_1 t$','interpreter','latex','fontsize',22)
xlim([-dt/10, td(end)+dt/10])
ylim([-L/300, L+L/300])


% figure()
% plot(td, sqrt(sum(error.^2, 2)) / input_dim, 'b-')

compute_lyap = 1;
num_lyap = 30;

[td, tt, xx, u_true, prediction, ly_next, ly_history] ...
            = TDRC_KS(m, set_average_degree, approx_reservoir_size,...
                      rho, num_delay, gama, beta, c, num_inter_step, ...
                      compute_lyap, num_lyap);
                  
figure('name', 'KS Lyapunov', 'position', [500,400,800,700])
plot(ly_history(:,1), 'linewidth', 2);hold on
plot(ly_history(:,2), 'linewidth', 2);hold on
plot(ly_history(:,3), 'linewidth', 2);hold on
plot(ly_history(:,4), 'linewidth', 2);hold on
plot(ly_history(:,5), 'linewidth', 2);hold on
set(gca,'Position',[0.1,0.11,0.85,0.85],'fontsize',20)
xlabel('Iterations', 'fontsize', 22);
ylabel('$\Lambda_j$', 'interpreter', 'latex', 'fontsize', 22)

save(['TDRC_KS_nonhomo', '_N', num2str(approx_reservoir_size), '_delay', num2str(num_delay), '.mat'],...
      'dt','approx_reservoir_size','num_delay','rho','gama','beta', 'c','num_inter_step', 'num_lyap', 'u_true','td', 'tt','xx', 'prediction', 'ly_history')





