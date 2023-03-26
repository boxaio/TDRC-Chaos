clc,clear;

% This program call function 'KS_Solver' to plot the solution of KS
% equation, and call function 'KS_Lyap' to compute the Lyapunov exponents of KS equation 

L = 60;
N = 128;
dt = 0.25;
tmax = 25000;

% mu = 0 for the standard spatially homogeneous KS equation
mu = 0.0;

% spatial nodes
x = L * (-N/2 + 1 : N/2)' / N;
u0 = 0.6 * (-1 + 2*rand(size(x)));

% 
% ut = KS_Solver(L, N, dt, tmax, mu, u0);
% 
% 
% save(['KS_L' num2str(L) '_mu' num2str(mu) '.mat'], 'ut', 'L', '-v7.3');
% 
% 
% % Plot contour
% figure('name', 'KS solution', 'position', [200,250,1600,550])
% td = dt * (1 : 400);
% xd = x + L/2;
% [xx, tt] = meshgrid(xd, td);
% pcolor(tt, xx, ut(1:400,:))
% shading interp
% colormap(slanCM('spectral', 128))
% % colormap(slanCM('haline', 128))
% cb = colorbar;
% % set(cb,'YTick',[-2.5, -1.0, 0, 1.0, 2.5], 'fontsize',20)
% % set(cb,'YTickLabel',{'-2.5', '-1.0', '0', '1.0', '2.5'}) 
% xlabel('$t$','interpreter','latex','fontsize',24)
% ylabel('$x$','interpreter','latex','fontsize',24)  
% set(gca,'Position',[0.06,0.15,0.88,0.8],'fontsize',24)
% set(gca,'fontsize',24)


mu = 0.1;
tmax = 4000; 
num_lyaps = 40; 
norm_steps = 10;

complex_spectrum = KS_Lyap(L, N, dt, tmax, mu, u0, num_lyaps, norm_steps);

spectrum = real(complex_spectrum);

spec_sum = zeros(num_lyaps, 1);

for i = 1 : num_lyaps
    spec_sum(i) = sum(spectrum(1:i));
end

abs_spec_sum = abs(spec_sum);
ky_point = min(abs_spec_sum);
kydim = find(spec_sum < 0, 1, 'first');


figure('name', 'KS Lyapunov', 'position', [500,400,800,700])
p1 = plot(spectrum(1:26), '+', 'linewidth', 2); hold on
set(get(p1, 'parent'), 'linewidth', 1.9)
p2 = refline([0,0]);
set(p2, 'color', 'k', 'linewidth', 2)
set(gca,'Position',[0.12,0.11,0.85,0.85],'fontsize',20)
xlabel('$j$', 'interpreter', 'latex','fontsize', 22);
ylabel('$\Lambda_j$', 'interpreter', 'latex', 'fontsize', 22)
