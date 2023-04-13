function complex_spectrum = KS_Lyap(L, N, dt, tmax, mu, u0, num_lyaps, norm_steps)
% Computes the Lyapunov Exponents of the Kuramoto-Sivashinsky equation with
% a multiplicative asymmetry term that breaks spatial homogeneity.

    
x = L * (-N/2+1 : N/2)' / N;

rng('shuffle')

% sets the spatial inhomgeneity wavelength
wavelength = L / 4; 
omega = 2 * pi / wavelength;
p = mu.*cos(omega.*x);
px = -omega*mu.*sin(omega.*x);
pxx = -(omega^2).*p;

% initial condition
v = fft(u0);

% tangent vectors
Y = orth(rand(N, num_lyaps));

% wave numbers
k = [0:N/2-1 0 -N/2+1:-1]' * (2 * pi / L);
% Fourier multipliers
L = k.^2 - k.^4;
E = exp(dt * L); E2 = exp(dt * L / 2);
M = 16;    % no. of points for complex means    
r = exp(1i * pi * ((1 : M) - 0.5) / M);   % roots of unity
LR = dt * L(:,ones(M,1)) + r(ones(N,1),:);

Q = dt * real(mean( (exp(LR/2) - 1)./LR ,2));
f1 = dt * real(mean( (-4 - LR + exp(LR).*(4 - 3 * LR + LR.^2))./LR.^3 ,2));
f2 = dt * real(mean( (2 + LR + exp(LR).*(-2 + LR))./LR.^3 ,2));
f3 = dt * real(mean( (-4 - 3 * LR - LR.^2 + exp(LR).*(4 - LR))./LR.^3 ,2));

nmax = round(tmax / dt); 
g = -0.5i * k;

vv = zeros(N, nmax);
Rii = zeros(num_lyaps, nmax / norm_steps);

vv(:,1) = v;
transient = 1000;
for n = 1 : transient
    r_ifft_v = real(ifft(v));
    Nv = g.*fft(r_ifft_v.^2) + 2i * k.*fft(r_ifft_v.*px) ...
          - fft(r_ifft_v.*pxx) + k.^2.*fft(r_ifft_v.*p);
    a = E2.*v + Q.*Nv;
    r_ifft_a = real(ifft(a));
    Na = g.*fft(r_ifft_a.^2) + 2i*k.*fft(r_ifft_a.*px) ...
          - fft(r_ifft_a.*pxx) + k.^2.*fft(r_ifft_a.*p);
    b = E2.*v + Q.*Na;
    r_ifft_b = real(ifft(b));
    Nb = g.*fft(r_ifft_b.^2) + 2i*k.*fft(r_ifft_b.*px) ...
          - fft(r_ifft_b.*pxx) + k.^2.*fft(r_ifft_b.*p);
    c = E2.*a + Q.*(2 * Nb - Nv);
    r_ifft_c = real(ifft(c));
    Nc = g.*fft(r_ifft_c.^2) + 2i.*k.*fft(r_ifft_c.*px) ...
          - fft(r_ifft_c.*pxx) + k.^2.*fft(r_ifft_c.*p);
    v = E.*v + Nv.*f1 + 2 * (Na + Nb).*f2 + Nc.*f3;
end

count = 0;


for n = 1 : nmax

    % Evolve KS 
   
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
    Nc = g.*fft(r_ifft_c.^2) + 2i.*k.*fft(r_ifft_c.*px) ...
          - fft(r_ifft_c.*pxx) + k.^2.*fft(r_ifft_c.*p);
    v = E.*v + Nv.*f1 + 2 * (Na + Nb).*f2 + Nc.*f3;
    vv(:,n) = v;

    % Evolve tangent map of KS
    Hv = (-1i*k).*fft(real(ifft(v)).*real(ifft(Y))) ...
          + (1i*k).*(fft(px.*real(ifft(Y)))) ...
          - fft(pxx.*real(ifft(Y)));
    a_ = E2.*Y + Q.*Hv;
    Ha = (-1i*k).*fft(real(ifft(v)).*real(ifft(a_))) ...
          + (1i*k).*fft(px.*real(ifft(Y))) ...
          - fft(pxx.*real(ifft(Y)));
    b_ = E2.*Y + Q.*Ha;
    Hb = (-1i*k).*fft(real(ifft(v)).*real(ifft(b_))) ...
          + (1i*k).*fft(px.*real(ifft(Y))) ...
          - fft(pxx.*real(ifft(Y)));
    c_ = E2.*a_ + Q.*(2 * Hb - Hv);
    Hc = (-1i*k).*fft(real(ifft(v)).*real(ifft(c_))) ...
          + (1i*k).*fft(px.*real(ifft(Y))) ...
          - fft(pxx.*real(ifft(Y)));
    Y = E.*Y + Hv.*f1 + 2 * (Ha + Hb).*f2 + Hc.*f3;

    % Normalize tangent vectors and record normalization

    if mod(n, norm_steps) == 0
        count = count + 1;
        % QR decomposition, Q is of size(N, num_lyaps), R is of size(num_lyaps, num_lyaps)
        [matQ, matR] = qr(real(ifft(Y)), 0);
        Rii(:, count) = log(diag(matR));
        % new Y
        Y = fft(matQ(:, 1:num_lyaps));
    end

end

complex_spectrum = sum(Rii, 2)./(nmax * dt); 

% spectrum = real(complex_spectrum);
% 
% spec_sum = zeros(num_lyaps,1);
% 
% for i = 1 : num_lyaps
%     spec_sum(i) = sum(spectrum(1:i));
% end
% 
% abs_spec_sum = abs(spec_sum);
% ky_point = min(abs_spec_sum);
% kydim = find(spec_sum < 0, 1,'first');
% 
% 
% figure()
% plot(spectrum, '+')
% refline([0,0])
% filename = strcat('L', num2str(L),'asym_m_delta', num2str(100*mu), 'wl', num2str(wavelength), 'numl', num2str(num_lyaps));
% %save(filename, 'spectrum', 'spec_sum', 'ky_point', 'kydim', 'Rii', 'delta', 'wavelength');

