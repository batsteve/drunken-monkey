 load('data/data_4_re_2.dat', 'par', 'TT', 'e_diss_mu', 'f_coeff',...
     'e_inflow_mu', '-mat');
 
par;
figure(1);
clf;
%xx = 1024+(1:2^11);
xx= 1:(5*4990);
subplot(2, 1, 1)
plot(TT(xx), e_diss_mu(xx), 'LineWidth', 2);
title('dissipation')
ylabel('A')
subplot(2, 1, 2)
plot(TT(xx), -abs(f_coeff_1(xx, 8, 7)), 'LineWidth', 2);
title('(1, 0) Fourier mode')
ylabel('B')
xlabel('time')


filename = sprintf('pics/time_series.jpg');
print(sprintf('-f%d', 1), '-painters', filename,'-djpeg');
filename = sprintf('pics/time_series.fig');
savefig(figure(1), filename, 'compact')
