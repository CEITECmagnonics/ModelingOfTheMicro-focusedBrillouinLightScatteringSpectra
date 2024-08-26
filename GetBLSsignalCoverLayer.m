function [BLSsignal, ff] = GetBLSsignalCoverLayer(hNiFe, Bext, n, mu, hPt)
% clearvars;
% clc;
% clf;
% Used constants
c=2.998e8; %m s^-1

%Experimental constants
lambda=532e-9;
f=c./lambda;
k0=2*pi/lambda;
NA=0.75;
% hNiFe=26e-9;

FocalLength=1.65e-3;
BeamWaist=500e-9;


DFNiFe= -8.1653 + 1j*15.348;
DFPt = -23.492 + 1j*4.6848;
% DFNiFe = (2.2 + 1j*10e-6)^2; %YIG
ExtinCoefNiFe = sqrt((abs(DFNiFe)-real(DFNiFe))/2);
DFSi=17.237 + 1j*0.43004;
% DFSi = (1.9 + 1j*10e-6)^2; %GGG

DF=[1,DFPt,DFNiFe,DFSi];
PM=ones(size(DF));

N=60;
qxHalf=linspace(0,1.1,N)*k0;
qx=cat(2,-flip(qxHalf(2:end)),qxHalf);
dqx=qx(2:end)-qx(1:end-1);
dqx=cat(1,0,dqx(:),0);
dqx=(dqx(2:end)+dqx(1:end-1))/2;
% Nq=length(qx)^2;

qy=qx;
dqy=dqx;
% dQ=dqx*dqy';

[Qx,Qy]=ndgrid(qx,qy);
Q=sqrt(Qx.^2+Qy.^2);
% Phi=atan2(Qy,Qx);
% ks=sqrt(DFNiFe)*k0;
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
% fftEI(:,:,1) = interp2(KX, KY, fftE(:,:,1), Qx, Qy, 'linear', 0);
% fftEI(:,:,2) = interp2(KX, KY, fftE(:,:,2), Qx, Qy, 'linear', 0);
% fftEI(:,:,3) = interp2(KX, KY, fftE(:,:,3), Qx, Qy, 'linear', 0);

[~, ~, Exi, Eyi, Ezi] = focalFields(lambda, NA, 0, 2, FocalLength, pi/dqx(2), 2*N-1);

fftEI(:,:,2) = fftshift(fftshift(fft2(ifftshift(ifftshift(Exi,1),2)),1),2);
fftEI(:,:,1) = fftshift(fftshift(fft2(ifftshift(ifftshift(Eyi,1),2)),1),2);
fftEI(:,:,3) = fftshift(fftshift(fft2(ifftshift(ifftshift(Ezi,1),2)),1),2);

zs = linspace(0, hNiFe, 100);
Volume = exp(-1*ExtinCoefNiFe*k0*zs);
VolumeFac = trapz(zs, Volume);

fftEI = fftEI.*VolumeFac;

% fftEI(:,:,1) = fftEI(:,:,1)./max(fftEI(:));
% fftEI(:,:,2) = fftEI(:,:,2)./max(fftEI(:));
% fftEI(:,:,3) = fftEI(:,:,3)./max(fftEI(:));

[ff, f00fkx] = SpinWaveGreen(qx+0.0001, Bext, hNiFe, length(qx), 90, n, mu);
f00fkx(isnan(f00fkx)) = 0;

% load("M:\MICROMAG\Wojewoda\GreenFunctionNiFe\!Analysis\50mT-All3Comp-Complex\CorrectFFT.mat");
% load("M:\MICROMAG\Wojewoda\GreenFunctionNiFe\!Analysis\50mT-All3Comp-Complex\FFTz.mat");
% ff=fS;

ExS2 = zeros(size(ff));
BLSsignal = zeros(size(ff));
ExS = zeros(size(ff));
EyS = zeros(size(ff));
i = 0;
for fi=ff
    i=i+1;
    
    %Calculation of polarization   
    Mz = squeeze(f00fkx(i,:,:));
    My = 1j*squeeze(f00fkx(i,:,:));
    Mx = zeros(size(My));

%     FFTtruncYI = mkY(:,:,i);
%     FFTtruncZI = mkZ(:,:,i);
%     
%     [KXorig, KYorig] = meshgrid(kS, kS);
%     My = interp2(KXorig*1e6, KYorig*1e6, FFTtruncYI, Qx, Qy, 'linear', 0);
%     Mz = interp2(KXorig*1e6, KYorig*1e6, FFTtruncZI, Qx, Qy, 'linear', 0);
%     Mx = interp2(KXorig*1e6, KYorig*1e6, FFTtruncXI, Qx, Qy, 'linear', 0);
%     Mx = zeros(size(Mz));
    
%     Mx = Mx./max(max(Mx));
%     My = My./max(max(My));
%     Mz = Mz./max(max(Mz));

    Px = conv2(fftEI(:,:,3),1j*My) + conv2(fftEI(:,:,2),-1j*Mz);
    Py = conv2(fftEI(:,:,1),1j*Mz) + conv2(fftEI(:,:,3),-1j*Mx);
    Pz = conv2(fftEI(:,:,2),1j*Mx) + conv2(fftEI(:,:,1),-1j*My);

%     if i==34
%         pause
%     end
    Kx = linspace(2*min(qx), 2*max(qx), length(Px));
    Ky = Kx;

    %Interpolation of polarization
    fftPx=interpn(Kx,Ky,Px,Qx,Qy,'linear');
    fftPy=interpn(Kx,Ky,Py,Qx,Qy,'linear');
    fftPz=interpn(Kx,Ky,Pz,Qx,Qy,'linear');


    [htp,hts]=Fresnelq(lambda,DF(:),PM(:),[hPt, hNiFe],2,1);
    tp=htp(Q);
    ts=hts(Q);
    tp{1}(isnan(tp{1}))=0;
    ts{1}(isnan(ts{1}))=0;
    tp{2}(isnan(tp{2}))=0;
    ts{2}(isnan(ts{2}))=0;
    [pGF,sGF]=sphGreenFunction(Qx,Qy,DFNiFe,lambda,tp,ts);
    
    Ep=pGF{1,1}.*fftPx.*exp(-1i*Kzs.*hNiFe)+pGF{1,2}.*fftPx.*exp(+1i*Kzs.*hNiFe);
    Ep=Ep+pGF{2,1}.*fftPy.*exp(-1i*Kzs*hNiFe)+pGF{2,2}.*fftPy.*exp(+1i*Kzs*hNiFe);
    Ep=Ep+pGF{3,1}.*fftPz.*exp(-1i*Kzs*hNiFe)+pGF{3,2}.*fftPz.*exp(+1i*Kzs*hNiFe);

    Es=sGF{1,1}.*fftPx.*exp(-1i*Kzs*hNiFe)+sGF{1,2}.*fftPx.*exp(+1i*Kzs*hNiFe);
    Es=Es+sGF{2,1}.*fftPy.*exp(-1i*Kzs*hNiFe)+sGF{2,2}.*fftPy.*exp(+1i*Kzs*hNiFe);
    Es=Es+sGF{3,1}.*fftPz.*exp(-1i*Kzs*hNiFe)+sGF{3,2}.*fftPz.*exp(+1i*Kzs*hNiFe);
    
    
    cosPhi=Qx./Q;
    sinPhi=Qy./Q;
    cosPhi(isnan(cosPhi))=1;
    sinPhi(isnan(sinPhi))=0;

    Ex=VolumeFac*Ep.*cosPhi-VolumeFac*Es.*sinPhi;
    Ey=VolumeFac*Ep.*sinPhi+VolumeFac*Es.*cosPhi;

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

    Ex=(FocalLength/k0)^2*numel(Ex)*1/(4*pi^2)*dkx*dky*fftshift(fftshift(ifft2(ifftshift(ifftshift(Ex.*double(Q<=k0*NA),1),2)),1),2);
    Ey=(FocalLength/k0)^2*numel(Ey)*1/(4*pi^2)*dkx*dky*fftshift(fftshift(ifft2(ifftshift(ifftshift(Ey.*double(Q<=k0*NA),1),2)),1),2);

    Ex=Ex.*PSFFilter;
    Ey=Ey.*PSFFilter;
    

    ExS2(i)=dxi*dyi*sum(Ex(:).*conj(Ex(:)));
    BLSsignal(i)=dxi*dyi*sum(Ex(:).*conj(Ex(:)));
    
    ExS(i)=dxi*dyi*sum(Ex(:));
    EyS(i)=dxi*dyi*sum(Ey(:));
end

end






