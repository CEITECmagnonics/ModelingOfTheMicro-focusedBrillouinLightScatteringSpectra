hNiFe = 31.5E-9;
% Bexts = linspace(1e-3, 550e-3, 51);
Bexts = [10e-3, 558e-3];

i=0;
sumTime=0;
BLSsignaln0 = zeros(30,length(Bexts));
ffn0 = zeros(30,length(Bexts));
BLSsignaln1 = zeros(30,length(Bexts));
ffn1 = zeros(30,length(Bexts));
for Bext =Bexts
   i=i+1;
   timer = tic;   
   [BLSsignaln0(:,i), ffn0(:,i)] = GetBLSsignal(hNiFe, Bext, 0, -1000e9);
   [BLSsignaln1(:,i), ffn1(:,i)] = GetBLSsignal(hNiFe, Bext, 1, -1000e9);
   elTime = toc(timer);
   sumTime = sumTime + elTime;
   elAvgTime = sumTime/i;
   if mod(i,2)==0
    RemainTime = (length(Bexts)-i)*elAvgTime
   end
end



BextsMatrix = ones(size(ffn0)).*Bexts;



ffE = linspace(0, max(ffn1(:)), 201);
[BEXTi, ffEi] = ndgrid(Bexts, ffE);

i=0;
BLSsignalIn0 = zeros(size(BEXTi));
BLSsignalIn1 = zeros(size(BEXTi));
for Bext =Bexts
    i=i+1;
    BLSsignalIn0(i,:) = interp1(ffn0(:,i), BLSsignaln0(:,i), ffE, 'linear', 0);
    BLSsignalIn1(i,:) = interp1(ffn1(:,i), BLSsignaln1(:,i), ffE, 'linear', 0);
end

figure()
surf(BEXTi, ffEi, (BLSsignalIn0 + BLSsignalIn1), 'EdgeColor', 'None')
view(2)
axis tight

figure()
plot(squeeze(ffEi(1,:)), squeeze(BLSsignalIn0(1,:)))
hold on
plot(ffn0(:,1), squeeze(BLSsignaln0(:,1)))

figure()
plot(squeeze(ffEi(1,:)), squeeze(BLSsignalIn0(15,:)))
hold on
plot(ffn0(:,15), squeeze(BLSsignaln0(:,15)))


