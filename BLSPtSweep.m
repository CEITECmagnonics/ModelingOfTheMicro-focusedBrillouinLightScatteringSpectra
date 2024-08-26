hPts = linspace(0, 50e-9, 51);
i = 0;
sumTime=0;
multiWaitbar('main', 'Reset');
for hPti = hPts
   i=i+1;
%    timer = tic;   
   [BLSsignaln0(:,i), ffn0(:,i)] = GetBLSsignalCoverLayer(31.5e-9, 10e-3, 0, -1000e9,hPti);
%    elTime = toc(timer);
%    sumTime = sumTime + elTime;
%    elAvgTime = sumTime/i;
%    if mod(i,2)==0
%         RemainTime = (length(hPts)-i)*elAvgTime
%    end
       multiWaitbar('main', (i-1)/length(hPts));
end
    
[X, Y] = ndgrid(ffn0(:,i), hPts*1e9);
figure('name', 'BLS vs thickness of Pt')
surf(X*1e6, Y, BLSsignaln0, 'EdgeColor', 'None');
view(2)
xlabel('Frequency (GHz)')
ylabel('Thickness of Pt (nm)')