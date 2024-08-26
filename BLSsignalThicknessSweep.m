Bext = 500e-3;
hNiFes = linspace(1e-9, 100e-9, 101);

i=0;
sumTime=0;
BLSsignaln0 = zeros(30,length(hNiFes));
ffn0 = zeros(30,length(hNiFes));
BLSsignaln1 = zeros(30,length(hNiFes));
ffn1 = zeros(30,length(hNiFes));
multiWaitbar('main', 'Reset');
for hNiFe =hNiFes
   i=i+1;
%    timer = tic;   
   [BLSsignaln0(:,i), ffn0(:,i)] = GetBLSsignal(hNiFe, Bext, 0, -1000e9);
   [BLSsignaln1(:,i), ffn1(:,i)] = GetBLSsignal(hNiFe, Bext, 1, -1000e9);
   multiWaitbar('main', (i-1)/length(hNiFes));
%    elTime = toc(timer);
%    sumTime = sumTime + elTime;
%    elAvgTime = sumTime/i;
%    if mod(i,4)==0
%     RemainTime = (length(hNiFes)-i)*elAvgTime/60;
%     sprintf("Remaining time: %d minutes", round(RemainTime))
%    end
end



BextsMatrix = ones(size(ffn0)).*hNiFes;



ffE = linspace(20, 35, 201);
[hNIFEi, ffEi] = ndgrid(hNiFes, ffE);

i=0;
BLSsignalIn0 = zeros(size(hNIFEi));
BLSsignalIn1 = zeros(size(hNIFEi));
for hNiFe =hNiFes
    i=i+1;
    BLSsignalIn0(i,:) = interp1(ffn0(:,i), BLSsignaln0(:,i), ffE, 'linear', 0);
    BLSsignalIn1(i,:) = interp1(ffn1(:,i), BLSsignaln1(:,i), ffE, 'linear', 0);
end

figure()
surf(hNIFEi, ffEi, (BLSsignalIn0 + BLSsignalIn1), 'EdgeColor', 'None')
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


