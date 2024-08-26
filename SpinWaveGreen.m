% Created on Mon Oct 26
% This script calculates density of states for BLS signal calculations
% If you do not now how to call Python from Matlab start with 
% https://www.mathworks.com/help/matlab/matlab_external/create-object-from-python-class.html
% @author: Ondrej Wojewoda
function [ff, f00fkx] = SpinWaveGreen(kx, Bext, d, N, Nf, n, mu)
pathToSWT = fileparts(which('SpinWaveToolkit.py'));
if count(py.sys.path,pathToSWT) == 8
    insert(py.sys.path,int32(0),pathToSWT);
end


% kx = linspace(-25.001e6, 25e6, 150);
kx = kx*sqrt(2);
kxi = py.numpy.linspace(min(kx), max(kx), int32(length(kx)));

% Dispersion of two first thickness modes of 30 nm NiFe thin film
material = py.SpinWaveToolkit.Material(800e3, 16e-12, 70e-4,  0, 28.8*2*pi*1e9, 0);
% material = py.SpinWaveToolkit.Material(99e3, 4.2e-12, 13e-4,  0, 28*2*pi*1e9, 0);
% material = py.SpinWaveToolkit.YIG();
theta = pi/2;
phii = linspace(0, pi, 200);
weff = 3e-6;
boundaryCond = 1;
% Bext = 0.55;
% d = 30e-9;

f00P = zeros(length(phii), length(kx));
lifetimeP = zeros(length(phii), length(kx));

% f00PSAFM = zeros(length(phii), 4, length(kx));
% lifetimeSAFM = zeros(length(phii), 4, length(kx));
% lifetimeP = zeros(length(phii), 4, length(kx));

i=0;
for phi=phii
    i=i+1;
    NiFeChar = py.SpinWaveToolkit.DispersionCharacteristic(Bext, material, d, kxi, theta, phi, weff, boundaryCond);
    f00P(i,:) =  double(NiFeChar.GetDispersion(n))*1e-9/(2*pi);
    lifetimeP(i,:) =  double(NiFeChar.GetLifetime(n))*1e9;

%     CoFeB = py.SpinWaveToolkit.DispersionCharacteristic(0, py.SpinWaveToolkit.CoFeB(), 15e-9, kxi, theta, phi, 0, boundaryCond,0, 1.5e3, 1.5e3, -0.667e-3, -0.284e-3, 0.6e-9, 15e-9, py.SpinWaveToolkit.CoFeB(), -0.667e-3, -0.284e-3, phi);
%     f00PSAFM(i,:,:) =  double(CoFeB.GetDispersionSAFMNumeric())*1e-9/(2*pi);
%     lifetimeP(i,:)  = 2*1e9;
%     f00P(i,:) = f00PSAFM(i,3,:); 
end
f00P = squeeze(f00P);
[RHO, PHI] = meshgrid(kx, phii);
[X, Y] = pol2cart(PHI, RHO);
xi = linspace(-max(kx)/sqrt(2), max(kx)/sqrt(2), N);
yi = linspace(-max(kx)/sqrt(2), max(kx)/sqrt(2), N);
[Xi, Yi] = ndgrid(xi, yi);

Ff00 = scatteredInterpolant(X(:), Y(:), f00P(:), 'linear', 'none');
f00 = Ff00(Xi, Yi);
Flifetime = scatteredInterpolant(X(:), Y(:), lifetimeP(:), 'linear', 'none');
lifetime = Flifetime(Xi, Yi);

% figure('name', 'Dispersion relation')
% plot(kxi*1e-6, f00, kxi*1e-6, f11);
% xlabel('kxi (rad/um)');
% ylabel('Frequency (GHz)');
% legend('00', '11')
% title('Dispersion relation of NiFe n=0,1')


ff = linspace(min(f00(:))-2, max(f00(:))+2, Nf);
f00fkx = zeros(length(ff),N,N);
i=0;
% figure()


% [Kx, Ky] = ndgrid(xi/1e6, yi/1e6);
for f=ff
    i=i+1;
    w = 2*pi*f;
    w00 = 2*pi*f00;
    % Bose-Einstein distribution
    hbar = 1.0545718e-34;
    kb = 1.38064852e-23;
    T = 300;
    mu=mu;
    B=1;
    BE = 1./(B*exp((hbar*2*pi*(abs(f*1e9)-mu))/(kb*T))-1);
    RJ = hbar*2*pi*(abs(f*1e9)./(hbar*2*pi*(abs(f*1e9)-mu)));
%     f00fkx(i,:,:) = ((BE.^2*(w00.^2+w.^2))./(abs(w00.^2-w.^2).^2+(4./lifetime).^2.*w.^2) + (BE.^2*(w00.^2+w.^2))./(abs(w00.^2+w.^2).^2+(4./lifetime).^2.*w.^2));
    f00fkx(i,:,:) = BE./(abs(w00-w).^2+(2./lifetime).^2);
%     f00fkx(i,:,:) = BE./(4*pi^2*abs(f00-f).^2+(2./lifetime).^2); 
        
%     sgtitle(sprintf('Frequency: %.2f GHz',f))
%     ax = subplot(1,2,1);  
%     hold off
%     surf(Xi./1e6, Yi./1e6, f00, 'EdgeColor', 'None')
%     hold on
%     surf(Xi./1e6, Yi./1e6, f*ones(size(f00)), 'FaceAlpha', 0.8, 'EdgeColor', 'None') 
%     xlabel('k_x (\mum)')
%     ylabel('k_y (\mum)')
%     zlabel('Frequency (GHz)')
%     ax.DataAspectRatio = [1 1 0.12];
%     subplot(1,2,2)
%     surf(Xi./1e6, Yi./1e6, squeeze(f00fkx(i,:,:)), 'EdgeColor', 'None')
%     xlabel('k_x (\mum)')
%     ylabel('k_y (\mum)')
%     grid off
%     axis equal
%     fg = gcf;
%     if i==90
%         exportgraphics(fg,sprintf('DispGif%d.png', i),'Resolution',300)
%     end
end
end

