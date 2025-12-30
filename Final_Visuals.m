% FINAL_VISUALS

figure(10); clf;
set(gcf, 'Color', 'w', 'Position', [100 100 1200 800]);

% ---density--
rho = real(U(:,:,1));

subplot(2,1,1);
pcolor(METRICS.Xc, METRICS.Yc, rho);
shading interp; colormap(jet); colorbar;

caxis([0.0, 0.5]); 

ylim([0, 0.8]); 
xlim([0, 2.7]);

title('Density (Zoomed on Engine Core)'); 
xlabel('Axial (m)'); ylabel('Radial (m)');

% --- mach ---
rho_safe = rho + 1e-16;
u = real(U(:,:,2)) ./ rho_safe;
v = real(U(:,:,3)) ./ rho_safe;
p = abs(real((phys.gamma - 1) * (U(:,:,4) - 0.5 * rho .* (u.^2 + v.^2))));
a = sqrt(phys.gamma * p ./ rho_safe);
Mach = real(sqrt(u.^2 + v.^2) ./ a);

subplot(2,1,2);
pcolor(METRICS.Xc, METRICS.Yc, Mach);
shading interp; colormap(jet); colorbar;

caxis([1.5, 4.0]); 

ylim([0, 0.8]);
xlim([0, 2.7]);

title('Mach Number (Zoomed on Engine Core)'); 
xlabel('Axial (m)'); ylabel('Radial (m)');