function ut = KS_Solver(L, N, dt, tmax, mu, u0)
% % This program solve the Kuramoto-Sivashinsky equation:
% %    u_t + para(1) * u_xx + para(2) * u_xxxx + para(3) * u * u_x = 0
% % with u(x,0) = f(x),  u(0,t) = u(L,t)
% % using ode solver after a spatial discretization.

% Kassam A K, Trefethen L N. 
% Fourth-order time-stepping for stiff PDEs[J]. 
% SIAM Journal on Scientific Computing, 2005, 26(4): 1214-1233.

% spatial nodes
x = L*(-N/2+1:N/2)'/N;

wavelength = L / 4;
omega = 2 * pi / wavelength;
p = mu.*cos(omega.*x);
px = -omega*mu.*sin(omega.*x);
pxx = -(omega^2).*p;

% initial condition
v = fft(u0);

% wave numbers
k = [0:N/2-1 0 -N/2+1:-1]' * (2 * pi / L); 
% Fourier multipliers
L = k.^2 - k.^4;
E = exp(dt*L); E2 = exp(dt*L/2);
M = 16;   % no. of points for complex means
r = exp(1i * pi * ((1 : M) - 0.5) / M);  % roots of unity

LR = dt * L(:,ones(M, 1)) + r(ones(N, 1),:);

Q = dt * real(mean( (exp(LR/2)-1)./LR ,2));
f1 = dt * real(mean( (-4 - LR + exp(LR).*(4 - 3 * LR + LR.^2))./LR.^3, 2));
f2 = dt * real(mean( (2 + LR + exp(LR).*(-2 + LR))./LR.^3, 2));
f3 = dt * real(mean( (-4 - 3 * LR - LR.^2 + exp(LR).*(4 - LR))./LR.^3, 2));

nmax = round(tmax / dt); 
g = -0.5i * k;

vv = zeros(N, nmax);

vv(:, 1) = v;

for n = 1 : nmax
    r_ifft_v = real(ifft(v));
    Nv = g.*fft(r_ifft_v.^2) + 2i * k.*fft(r_ifft_v.*px) ...
          - fft(r_ifft_v.*pxx) + k.^2.*fft(r_ifft_v.*p);
    a = E2.*v + Q.*Nv;
    r_ifft_a = real(ifft(a));
    Na = g.*fft(r_ifft_a.^2) + 2i * k.*fft(r_ifft_a.*px) ...
          - fft(r_ifft_a.*pxx) + k.^2.*fft(r_ifft_a.*p);
    b = E2.*v + Q.*Na;
    r_ifft_b = real(ifft(b));
    Nb = g.*fft(r_ifft_b.^2) + 2i * k.*fft(r_ifft_b.*px) ...
          - fft(r_ifft_b.*pxx) + k.^2.*fft(r_ifft_b.*p);
    c = E2.*a + Q.*(2 * Nb - Nv);
    r_ifft_c = real(ifft(c));
    Nc = g.*fft(r_ifft_c.^2) + 2i * k.*fft(r_ifft_c.*px) ...
          - fft(r_ifft_c.*pxx) + k.^2.*fft(r_ifft_c.*p);
    v = E.*v + Nv.*f1 + 2 * (Na + Nb).*f2 + Nc.*f3;
    vv(:,n) = v;
end

ut = transpose(real(ifft(vv)));



