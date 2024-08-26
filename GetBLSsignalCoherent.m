function [ExS, EyS] =  GetBLSsignalCoherent(kS, NA, fil, B2)
% clearvars;
% clc;
% clf;
% Used constants
c=2.998e8; %m s^-1

%Experimental constants
lambda=532e-9;
f=c./lambda;
w=2*pi*f;
k0=2*pi/lambda;
% NA=0.75;
hNiFe=30e-9;

FocalLength=1.65e-3;
BeamWaist=1500e-9;


DFNiFe= -8.1653 + 1j*15.348;
DFSi=17.237 + 1j*0.43004;
% DFSi = (1.9 + 1j*10e-6)^2; %GGG
% DFNiFe = (2.2 + 1j*10e-6)^2; %YIG

DF=[1,DFNiFe,DFSi];
PM=ones(size(DF));

N=250;
qx=linspace(0,2,N)*k0;
qx=cat(2,-flip(qx(2:end)),qx);
dqx=qx(2:end)-qx(1:end-1);
dqx=cat(1,0,dqx(:),0);
dqx=(dqx(2:end)+dqx(1:end-1))/2;
Nq=length(qx)^2;

qy=qx;
dqy=dqx;
dQ=dqx*dqy';

[Qx,Qy]=ndgrid(qx,qy);
Q=sqrt(Qx.^2+Qy.^2);
Phi=atan2(Qy,Qx);
ks=sqrt(DFNiFe)*k0;
Kzs=(DFNiFe*k0^2-Q.^2).^(1/2);
Kz=(k0^2-Q.^2).^(1/2);

% This is for testing, and showing how the source which radiate to all
% direction behaves
% Px=zeros(size(Q));
% Py=ones(size(Q));
% Pz=ones(size(Q));
% Kx=Qx;
% Ky=Qy;

% [KX, KY] = meshgrid(kx, ky);
% fftEI = zeros(length(Qx), length(Qy), 3);
% fftEI(:,:,1) = interp2(KX, KY, fftE(:,:,2), Qx, Qy, 'linear', 0);
% fftEI(:,:,2) = interp2(KX, KY, fftE(:,:,1), Qx, Qy, 'linear', 0);
% fftEI(:,:,3) = interp2(KX, KY, fftE(:,:,3), Qx, Qy, 'linear', 0);

% [xi, yi, Exi, Eyi, Ezi] = focalFieldsRadial(lambda, NA, 0, fil, FocalLength, pi/dqx(2), 2*N-1);
[xi, yi, Exi, Eyi, Ezi] = focalFields(lambda, NA, 0, fil, FocalLength, pi/dqx(2), 2*N-1);
% Eyi = Eyi./sum(sum(sqrt(Exi.*conj(Exi) + Eyi.*conj(Eyi))));
% Exi = Exi./sum(sum(sqrt(Exi.*conj(Exi) + Eyi.*conj(Eyi))));
% sum(sum(sqrt(Exi.*conj(Exi) + Eyi.*conj(Eyi))));
fftEI(:,:,1) = fftshift(fftshift(fft2(ifftshift(ifftshift(Exi,1),2)),1),2);
fftEI(:,:,2) = fftshift(fftshift(fft2(ifftshift(ifftshift(Eyi,1),2)),1),2);
fftEI(:,:,3) = fftshift(fftshift(fft2(ifftshift(ifftshift(Ezi,1),2)),1),2);

% figure()
% surf(Qx/1e6, Qy/1e6, Exi.*conj(Exi), 'EdgeColor', 'None');
% view(2)

% kS=linspace(0, 20e6, 14);
[~, Ind0] = min(abs(qx-0));

i=1;
ExS = zeros(size(kS));
EyS = zeros(size(kS));
for k=kS
    [~, Ind] = min(abs(qx-k));
    Mx = zeros(size(Qx));
    My = zeros(size(Qx));
    Mz = zeros(size(Qx));
    
    Mx(Ind0, Ind) = 1;
    Mz(Ind0, Ind) = -1*1j;
    

    Px = conv2(fftEI(:,:,3),1j*My) + conv2(fftEI(:,:,2),-1j*Mz);
    Py = conv2(fftEI(:,:,1),1j*Mz) + conv2(fftEI(:,:,3),-1j*Mx);
    Pz = conv2(fftEI(:,:,2),1j*Mx) + conv2(fftEI(:,:,1),-1j*My);
    
%     Px = conv2(fftEI(:,:,3),1j*My) + conv2(fftEI(:,:,2),-1j*Mz) + B2*conv2(fftEI(:,:,3),Mz) + B2*conv2(fftEI(:,:,2),My);
%     Py = conv2(fftEI(:,:,1),1j*Mz) + conv2(fftEI(:,:,3),-1j*Mx) + B2*conv2(fftEI(:,:,1),My);
%     Pz = conv2(fftEI(:,:,2),1j*Mx) + conv2(fftEI(:,:,1),-1j*My) + B2*conv2(fftEI(:,:,1),Mz);


    Kx = linspace(2*min(qx), 2*max(qx), length(Px));

    Ky = Kx;

    %Interpolation of polarization
    fftPx=interpn(Kx,Ky,Px,Qx,Qy,'linear');
    fftPy=interpn(Kx,Ky,Py,Qx,Qy,'linear');
    fftPz=interpn(Kx,Ky,Pz,Qx,Qy,'linear');


    [htp,hts]=Fresnelq(lambda,DF(:),PM(:),[hNiFe],2,1);
    tp=htp(Q);
    ts=hts(Q);
    [pGF,sGF]=sphGreenFunction(Qx,Qy,DFNiFe,lambda,tp,ts);

    Ep=pGF{1,1}.*fftPx.*exp(-1i*Kzs*hNiFe)+pGF{1,2}.*fftPx.*exp(+1i*Kzs*hNiFe);
    Ep=Ep+pGF{2,1}.*fftPy.*exp(-1i*Kzs*hNiFe)+pGF{2,2}.*fftPy.*exp(+1i*Kzs*hNiFe);
    Ep=Ep+pGF{3,1}.*fftPz.*exp(-1i*Kzs*hNiFe)+pGF{3,2}.*fftPz.*exp(+1i*Kzs*hNiFe);

    Es=sGF{1,1}.*fftPx.*exp(-1i*Kzs*hNiFe)+sGF{1,2}.*fftPx.*exp(+1i*Kzs*hNiFe);
    Es=Es+sGF{2,1}.*fftPy.*exp(-1i*Kzs*hNiFe)+sGF{2,2}.*fftPy.*exp(+1i*Kzs*hNiFe);
    Es=Es+sGF{3,1}.*fftPz.*exp(-1i*Kzs*hNiFe)+sGF{3,2}.*fftPz.*exp(+1i*Kzs*hNiFe);

    cosPhi=Qx./Q;
    sinPhi=Qy./Q;
    cosPhi(isnan(cosPhi))=1;
    sinPhi(isnan(sinPhi))=0;

    Ex=Ep.*cosPhi-Es.*sinPhi;
    Ey=Ep.*sinPhi+Es.*cosPhi;

    Factor=(-2i*pi*sqrt(Kz*k0))*exp(1i*k0*FocalLength)/FocalLength;
    Ex=Ex.*Factor;
    Ey=Ey.*Factor;

    dkx=qx(2)-qx(1);
    dky=qy(2)-qy(1);
    Nxi=length(qx);
    Nyi=length(qy);

    DXi=2*pi/dkx;
    DYi=2*pi/dky;
    dxi=DXi/Nxi;
    dyi=DYi/Nyi;
    xi=linspace(-(Nxi-1)/2,(Nxi-1)/2,Nxi)*dxi;
    yi=linspace(-(Nyi-1)/2,(Nyi-1)/2,Nyi)*dyi;
    [Xi,Yi]=ndgrid(xi,yi);

    PSFFilter=exp(-(Xi.^2+Yi.^2)/BeamWaist^2);
%     PSFFilter=abs(Exi);
%     PSFFilter=ones(size(Xi));
    
    Ex=(FocalLength/k0)^2*numel(Ex)*1/(4*pi^2)*dkx*dky*fftshift(fftshift(ifft2(ifftshift(ifftshift(Ex.*double(Q<=k0*NA),1),2)),1),2);
    Ey=(FocalLength/k0)^2*numel(Ey)*1/(4*pi^2)*dkx*dky*fftshift(fftshift(ifft2(ifftshift(ifftshift(Ey.*double(Q<=k0*NA),1),2)),1),2);

    Ex=Ex.*PSFFilter;
    Ey=Ey.*PSFFilter;

    ExS(i)=dxi*dyi*sum(Ex(:));
    EyS(i)=dxi*dyi*sum(Ey(:));
    i=i+1;
end
% figure()
% surf(Xi*1e6, Yi*1e6, abs(Ex), 'EdgeColor', 'None');
% view(2)

% Exp = EyS.*conj(EyS);
% figure()
% hold on
% plot(kS/1e6, abs(EyS).^2)
% xlabel("Wavenumber (rad/um)")
% ylabel("BLS signal (arb. units)")
% 
% t = linspace(0,1,51);
% jj=1:length(kS);
% figure()
% for j=jj
% plot(ExS(j).*exp(1j*t*2*pi), EyS(j).*exp(1j*t*2*pi))
% hold on
% end
% axis equal
% 
% legendStrings = "kf = " + string(kS(jj)/1e6);
% legend(legendStrings)
end


% 
% plot(kS, abs(EyS).^2, kS, abs(EySrad).^2)
% legend('Linear polarization', 'Radial polarization')
% xlabel('Wavenumber (rad/m)')
% ylabel('BLS sensitivity (arb. units)')


