f = linspace(0.05, 2, 51);
i = 1;
ExiX = zeros(201,length(f));
for fi = f
    [xi, yi, Exi, Eyi, Ezi] = focalFields(532e-9, 0.75, 0, fi, 5e-3, 3e-6);
    ExiX(:,i) = Exi(:, round(end/2)).*conj(Exi(:, round(end/2)));
%     ExiX(:,i) = Exi(:, round(end/2));
    i = i + 1
end
    
[X, F] = ndgrid(xi, f);
figure('name', 'filling factor vs Ex')
surf(X*1e6, F, abs(ExiX), 'EdgeColor', 'None');
view(2)
xlabel('x (\mum)')
ylabel('Filling factor ()')