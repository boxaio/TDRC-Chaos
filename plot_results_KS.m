clc,clear;

m1 = matfile('TDRC_KS_N1000_delay3');
m2 = matfile('TDRC_KS_N800_delay5');
m3 = matfile('TDRC_KS_nonhomo_N1000_delay3');
m4 = matfile('TDRC_KS_nonhomo_N800_delay5');

dt = m1.dt;
L = 22;

figure('position', [500,400,600,500])
p1 = plot(m1.ly_history(end,:),'b+','linewidth',2);hold on;
set(get(p1,'parent'),'linewidth',1.9)
p2 = plot(m2.ly_history(end,:),'rx','linewidth',2);hold on
p3 = plot(m3.ly_history(end,:),'mo','linewidth',2);hold on
p4 = plot(m4.ly_history(end,:),'c*','linewidth',2);hold on
lgd = legend([p1,p2,p3,p4], {'$N=1000,\tau=3\Delta t,\mu=0$', ...
                              '$N=800,\tau=5\Delta t,\mu=0$', ...
                              '$N=1000,\tau=3\Delta t,\mu=0.07$',...
                              '$N=800,\tau=5\Delta t,\mu=0.07$'}, ...
                              'interpreter','latex','fontsize', 20,'location', 'southwest','box', 'off');  
lgd.ItemTokenSize = [60,20];  % set the legend length
p5 = plot([0, 31], [0,0], 'k--', 'linewidth', 1.8);hold off
set(get(get(p5, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
set(gca,'Position',[0.17,0.17,0.75,0.75],'fontsize',22)
xlabel('$j$', 'interpreter', 'latex', 'fontsize', 24);
ylabel('$\Lambda_j$', 'interpreter', 'latex', 'fontsize', 24)
xlim([0, 31])
ylim([-5, 2])


% Plot contour
error1 = m1.u_true - m1.prediction;

figure('name', 'KS solution', 'position', [200,250,900,600])
subplot(311)
p1 = pcolor(m1.tt, m1.xx, m1.u_true);
shading interp
colormap(slanCM('spectral', 128))
colorbar('Ticks', [-2,-1,0,1,2], 'position', [0.94, 0.2, 0.02, 0.6]);
set(get(p1, 'parent'), 'linewidth', 2)
title('True','fontsize',20)
set(gca,'XTicklabel',[]);
ylabel('$x$','interpreter','latex','fontsize',22)  
set(gca,'Position',[0.11,0.72,0.8,0.21],'fontsize',22)
xlim([-dt/10, m1.td(1,end)+dt/10])
ylim([-L/300, L+L/300])
subplot(312)
p1 = pcolor(m1.tt, m1.xx, m1.prediction);
shading interp
set(get(p1, 'parent'), 'linewidth', 2)
colormap(slanCM('spectral', 128))
title('Prediction','fontsize',20)
set(gca,'XTicklabel',[]);
ylabel('$x$','interpreter','latex','fontsize',22)  
set(gca,'Position',[0.11,0.43,0.8,0.21],'fontsize',22)
xlim([-dt/10, m1.td(1,end)+dt/10])
ylim([-L/300, L+L/300])
subplot(313)
p1 = pcolor(m1.tt, m1.xx, error1);
shading interp
set(get(p1, 'parent'), 'linewidth', 2)
colormap(slanCM('spectral', 128))
title('Error','fontsize',20)
ylabel('$x$','interpreter','latex','fontsize',22)  
set(gca,'Position',[0.11,0.13,0.8,0.21],'fontsize',22)
xlabel('$\Lambda_1 t$','interpreter','latex','fontsize',22)
xlim([-dt/10, m1.td(1,end)+dt/10])
ylim([-L/300, L+L/300])


% Plot contour
error4 = m4.u_true - m4.prediction;

figure('name', 'KS solution', 'position', [600,250,900,600])
subplot(311)
p1 = pcolor(m4.tt, m4.xx, m4.u_true);
shading interp
colormap(slanCM('spectral', 128))
colorbar('Ticks', [-2,-1,0,1,2], 'position', [0.94, 0.2, 0.02, 0.6]);
set(get(p1, 'parent'), 'linewidth', 2)
title('True','fontsize',20)
set(gca,'XTicklabel',[]);
ylabel('$x$','interpreter','latex','fontsize',22)  
set(gca,'Position',[0.11,0.72,0.8,0.21],'fontsize',22)
xlim([-dt/10, m4.td(1,end)+dt/10])
ylim([-L/300, L+L/300])
subplot(312)
p1 = pcolor(m4.tt, m4.xx, m4.prediction);
shading interp
set(get(p1, 'parent'), 'linewidth', 2)
colormap(slanCM('spectral', 128))
title('Prediction','fontsize',20)
set(gca,'XTicklabel',[]);
ylabel('$x$','interpreter','latex','fontsize',22)  
set(gca,'Position',[0.11,0.43,0.8,0.21],'fontsize',22)
xlim([-dt/10, m4.td(1,end)+dt/10])
ylim([-L/300, L+L/300])
subplot(313)
p1 = pcolor(m4.tt, m4.xx, error4);
shading interp
set(get(p1, 'parent'), 'linewidth', 2)
colormap(slanCM('spectral', 128))
title('Error','fontsize',20)
ylabel('$x$','interpreter','latex','fontsize',22)  
set(gca,'Position',[0.11,0.13,0.8,0.21],'fontsize',22)
xlabel('$\Lambda_1 t$','interpreter','latex','fontsize',22)
xticks([0,2,4,6,8,10,12,14,16])
xticklabels({'0','2','4','6','8','10','12','14','16'})
xlim([-dt/10, m4.td(1,end)+dt/10])
ylim([-L/300, L+L/300])





















