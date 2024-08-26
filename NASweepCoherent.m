i=0;
sumTime=0;


NAs = linspace(0.1, 1, 51);
% NAs = 0.1;
kS = linspace(0, 20e6, 51);
fil = 2;

ExS = zeros(length(kS),length(NAs));
EyS = zeros(length(kS),length(NAs));

for NA =NAs
   i=i+1;
   timer = tic;   
   [ExS(:,i), EyS(:,i)] = GetBLSsignalCoherent(kS, NA, fil, 0);
   elTime = toc(timer);
   sumTime = sumTime + elTime;
   elAvgTime = sumTime/i;
   if mod(i,4)==0
    RemainTime = (length(NAs)-i)*elAvgTime/60;
    sprintf("Remaining time: %d minutes", round(RemainTime))
   end
end

figure();
Exp = EyS.*conj(EyS);
plot(kS./1e6, Exp)


i=1;
for NA = NAs
    if Exp(end,i)<max(Exp(:,i)./100)
        row(i) = find(Exp(:,i)<max(Exp(:,i)./100), 1);
    else
        row(i) = 51;
    end
    i=i+1;
end

rowExp = (kS(row)/1e6)';
figure()
plot(NAs, kS(row)/1e6)