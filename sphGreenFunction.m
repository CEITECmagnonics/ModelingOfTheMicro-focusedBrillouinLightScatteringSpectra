function [pGF,sGF]=sphGreenFunction(Kx,Ky,DFNiFe,lambda,tp,ts)
c=3e9;
mu0=4*pi*1e-7;

k0=2*pi/lambda;
w=c*k0;
ks=k0*sqrt(DFNiFe);
pGF=cell(3,2);
sGF=cell(3,2);

Kr=sqrt(Kx.^2+Ky.^2);
cosPhi=Kx./Kr;
cosPhi(isnan(cosPhi))=1;
sinPhi=Ky./Kr;
sinPhi(isnan(sinPhi))=0;

Kzs=(DFNiFe*k0^2-Kr.^2).^(1/2)+eps;

sGF{1,1}=-1i*w^2*mu0/2*sinPhi.*ts{1}./Kzs;
sGF{2,1}=1i*w^2*mu0/2*cosPhi.*ts{1}./Kzs;
sGF{3,1}=zeros(size(Kr));

sGF{1,2}=-1i*w^2*mu0/2*sinPhi.*ts{2}./Kzs;
sGF{2,2}=1i*w^2*mu0/2*cosPhi.*ts{2}./Kzs;
sGF{3,2}=zeros(size(Kr));

pGF{1,1}=1i*w^2*mu0/2*cosPhi.*tp{1}/ks;
pGF{2,1}=1i*w^2*mu0/2*sinPhi.*tp{1}/ks;
pGF{3,1}=-1i*w^2*mu0/2*tp{1}.*Kr./(Kzs*ks);

pGF{1,2}=-1i*w^2*mu0/2*cosPhi.*tp{2}/ks;
pGF{2,2}=-1i*w^2*mu0/2*sinPhi.*tp{2}/ks;
pGF{3,2}=-1i*w^2*mu0/2*tp{2}.*Kr./(Kzs*ks);








