clc,clear;

m1 = matfile('TDRC_KS_N1000_delay3');
m2 = matfile('TDRC_KS_N800_delay5');
% m3 = matfile('TDRC_KS_N700_delay5');

figure('position', [500,400,600,500])
p1 = plot(m1.ly_history(end,:),'b+','linewidth',2);hold on;
set(get(p1,'parent'),'linewidth',1.9)
p2 = plot(m2.ly_history(end,:),'rx','linewidth',2);hold on
p3 = plot(m3.ly_history(end,:),'mo','linewidth',2);hold on
lgd = legend([p1,p2,p3], {'$N=1200,\tau=2\Delta t$', '$N=800,\tau=3\Delta t$', '$N=700,\tau=5\Delta t$'}, 'interpreter','latex','fontsize', 22, 'box', 'off');  
lgd.ItemTokenSize = [60,20];  % set the legend length
p4 = plot([0, 31], [0,0], 'k--', 'linewidth', 1.8);hold off
set(get(get(p4, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
set(gca,'Position',[0.17,0.17,0.75,0.75],'fontsize',22)
xlabel('$j$', 'interpreter', 'latex', 'fontsize', 24);
ylabel('$\Lambda_j$', 'interpreter', 'latex', 'fontsize', 24)
xlim([0, 21])
ylim([-17, 17])