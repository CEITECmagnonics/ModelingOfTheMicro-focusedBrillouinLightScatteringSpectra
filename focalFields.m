% This function calculates the focal fields of the given (objective) lens
% Lambda - wavelength of the used light
% NA - numeric aperture of the objective lens
% z - defocus of the beam (z=0 completely in focus)
% f0 - filling factor
% f - focal length of the objective lens
% N - number of the points in final field
% -------Outputs-------
% xi, yi - vectors to generate mesh by ndgrid
% Exi, Eyi, Ezi - Electric field complex vector in the desired plane -  E(x,y)
function [xi, yi, Exi, Eyi, Ezi] = focalFields(lambda, NA, z, f0, f, rhoMax, N)


E0 = 1; % Amplitude of the wave before focus

n1 = 1; %Index of refraction in first medium  - only changes magnitude
n2 = 1; %Index of refraction in second medium - only changes magnitude

k0=2*pi/lambda*n2;

thetaMax = asin(NA/n2);

% z = 0e-6; % 0 means focal plane
% w0 = 5e-3; % beam radius
% f = 10e-3; % focal length

theta = linspace(0, thetaMax, 41);

% f0 = w0/(f*sin(thetaMax)); % filling factor
fw = exp(-1/(f0^2)*(sin(theta).^2)/(sin(thetaMax)^2)); % apodization function

i=1;
phi = linspace(0, 2*pi, 45);
rho = linspace(1e-12, rhoMax, 180);
I00 = zeros(size(rho));
I01 = zeros(size(rho));
I02 = zeros(size(rho));

Ex = zeros(length(rho), length(phi));
Ey = zeros(length(rho), length(phi));
Ez = zeros(length(rho), length(phi));

for rhoi = rho
    I00(i) = trapz(fw.*cos(theta).^(1/2).*sin(theta).*(1+cos(theta)).*besselj(0,k0*rhoi*sin(theta)).*exp(1j*k0*z.*cos(theta)));
    I01(i) = trapz(fw.*cos(theta).^(1/2).*(sin(theta).^2).*besselj(1,k0*rhoi*sin(theta)).*exp(1j*k0*z.*cos(theta)));
    I02(i) = trapz(fw.*cos(theta).^(1/2).*sin(theta).*(1-cos(theta)).*besselj(2,k0*rhoi*sin(theta)).*exp(1j*k0*z.*cos(theta)));
    j=1;
    for phii = phi
        Ex(i,j) = 1j*k0*f/2*sqrt(n1/n2)*E0*exp(1j*k0*f)*(I00(i) + I02(i)*cos(2*phii));
        Ey(i,j) = 1j*k0*f/2*sqrt(n1/n2)*E0*exp(1j*k0*f)*I02(i)*sin(2*phii);
        Ez(i,j) = 1j*k0*f/2*sqrt(n1/n2)*E0*exp(1j*k0*f)*(-2*1j*I01(i)*sin(phii));
        j=j+1;
    end
    i=i+1;
end

[PHI, RHO] = meshgrid(phi, rho);
[X, Y] = pol2cart(PHI, RHO);
xi = linspace(min(min(X)), max(max(X)), N);
yi = linspace(min(min(Y)), max(max(Y)), N);
[Xi, Yi] = ndgrid(xi, yi);

% [XYu, ind] = unique([X(:), Y(:)], 'rows');

Fx = scatteredInterpolant(X(:), Y(:), Ex(:), 'linear', 'nearest');
Exi = Fx(Xi, Yi);
Fy = scatteredInterpolant(X(:), Y(:), Ey(:), 'linear', 'nearest');
Eyi = Fy(Xi, Yi);
Fz = scatteredInterpolant(X(:), Y(:), Ez(:), 'linear', 'nearest');
Ezi = Fz(Xi, Yi);
% 
% Plotting of the results
% figure('name', 'Beam spots')
% subplot(1,3,1)
% surf(Xi, Yi, real(Exi), 'EdgeColor', 'None')
% title('Ex')
% view(2)
% axis square
% subplot(1,3,2)
% surf(Xi, Yi, real(Eyi), 'EdgeColor', 'None')
% title('Ey')
% view(2)
% axis square
% subplot(1,3,3)
% surf(Xi, Yi, real(Ezi), 'EdgeColor', 'None')
% title('Ez')
% view(2)
% axis square
% pause
% % Plot the phase
% figure('name', 'Beam spots phase')
% subplot(1,3,1)
% surf(Xi, Yi, angle(Exi), 'EdgeColor', 'None')
% title('Ex')
% view(2)
% axis square
% subplot(1,3,2)
% surf(Xi, Yi, angle(Eyi), 'EdgeColor', 'None')
% title('Ey')
% view(2)
% axis square
% subplot(1,3,3)
% surf(Xi, Yi, angle(Ezi), 'EdgeColor', 'None')
% title('Ez')
% view(2)
% axis square
% % 
% % 
% % 
% figure('name', 'Cuts of Ex')
% plot(xi, real(squeeze(Exi(round(end/2),:))), xi, real(squeeze(Exi(:,round(end/2)))))

