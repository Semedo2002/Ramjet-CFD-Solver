% CHECK DATA RANGE
rho = real(U(:,:,1));
Mach = sqrt(U(:,:,2).^2 + U(:,:,3).^2) ./ sqrt(phys.gamma * (phys.gamma-1)*(U(:,:,4)-0.5*rho.*(U(:,:,2).^2+U(:,:,3).^2)./rho.^2)./rho);

fprintf('--- FIELD STATISTICS ---\n');
fprintf('Density:  Min = %.4f  | Max = %.4f\n', min(rho(:)), max(rho(:)));
fprintf('Mach:     Min = %.4f  | Max = %.4f\n', min(Mach(:)), max(Mach(:)));

if (max(rho(:)) - min(rho(:))) < 0.01
    fprintf('WARNING: The field is mathematically flat. The solver is not generating gradients.\n');
else
    fprintf('STATUS: Gradients exist. Adjust color limits to see them.\n');
end