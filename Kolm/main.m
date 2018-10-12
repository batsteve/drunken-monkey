clc;
clear;
global n1 n2 L1 L2 kx ky x y nu

par = Kol_Params();
par.t_start = 1e0;
par.t_end = 2e1;	% increase for longer data series
par.save_vector_fields = false;
par.plot_intermediate_vector_fields = false;
par.t_step = 1;
par.re = 40;		% increase for moar kaos!
par.nu = 1/par.re;
par.data_folder = 'data/test';
mkdir(par.data_folder)
par.indicator_filename = sprintf('%s/test.dat', par.data_folder);

L1=2*pi; n1=128;
L2=2*pi; n2=128;

[x,y,kx,ky] = gvars(par.n1,par.n2,par.L1,par.L2);
nu = par.nu;

[u1,u2] = u0(par);
[u1,u2] = dns2d(u1,u2,par);

