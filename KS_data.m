clc,clear;

% This program call function 'KS_Solver' to plot the solution of KS
% equation, and call function 'KS_Lyap' to compute the Lyapunov exponents of KS equation 

L = 22;
N = 64;
dt = 0.1;
tmax = 10000;

% mu = 0 for the standard spatially homogeneous KS equation
% mu = 0.0;
mu = 0.07;

% spatial nodes
x = L * (-N/2 + 1 : N/2)' / N;
u0 = 1 * (-1 + 2*rand(size(x)));

disp(['solving KS PDE...'])
% solution, size (length(t), N)
ut = KS_Solver(L, N, dt, tmax, mu, u0);


% Plot contour
figure('name', 'KS solution', 'position', [200,250,1200,550])
td = dt * (1 : 2000);
xd = x + L/2;
[xx, tt] = meshgrid(xd, td);
pcolor(tt, xx, ut(1:2000,:))
shading interp
colormap(slanCM('spectral', 128))
% colormap(slanCM('haline', 128))
cb = colorbar;
% set(cb,'YTick',[-2.5, -1.0, 0, 1.0, 2.5], 'fontsize',20)
% set(cb,'YTickLabel',{'-2.5', '-1.0', '0', '1.0', '2.5'}) 
xlabel('$t$','interpreter','latex','fontsize',24)
ylabel('$x$','interpreter','latex','fontsize',24)  
set(gca,'Position',[0.08,0.15,0.83,0.8],'fontsize',22)


num_lyaps = 30; 
norm_steps = 10;

disp(['calculating KS Lyapunov exponents...'])

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
p1 = plot(spectrum(1:16), '+', 'linewidth', 2); hold on
set(get(p1, 'parent'), 'linewidth', 1.9)
plot([0,20], [0,0], 'k-', 'linewidth', 1.8);
set(gca,'Position',[0.12,0.11,0.85,0.85],'fontsize',20)
xlabel('$j$', 'interpreter', 'latex','fontsize', 22);
ylabel('$\Lambda_j$', 'interpreter', 'latex', 'fontsize', 22)
% ylim([-0.4, 0.15])


disp(['save KS results for : '])
disp(['L = ',num2str(L)])
disp(['mu = ',num2str(mu)])
disp(['number of Lyapunov exponents : ', num2str(num_lyaps)])

% save(['KS_L' num2str(L) '_mu00', '.mat'], 'mu', 'ut', 'L', 'dt', 'spectrum');
save(['KS_L' num2str(L) '_mu007', '.mat'], 'mu', 'ut', 'L', 'dt', 'spectrum');
