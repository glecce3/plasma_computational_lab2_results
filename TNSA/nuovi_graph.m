
 n = [0.1 1 10]

T_ele = [3.7, 6.61, 2.21]
T_ion = [0.12, 0.75, 0.56]
E_ele = [0, 2.68, 17.51, 30.76, 22.52, 20.49]

figure
plot(n, T_ele, 'o--', 'LineWidth', 1.5)
hold on
plot(n, T_ion, 'o--', 'LineWidth', 1.5)
ylabel('Temperature [MeV]')
xlabel('Density [n_c]')
set(gca, 'XScale', 'log') 
title('Estimated electrons and ions temperatures at different plasma densities')
grid on
lgd = legend('T_{electrons}', 'T_{ions}') 
set(lgd, 'FontSize', 13);           % imposta la dimensione del font
