clear all; clc;

L = 16; v = 16/160;
n = 160;    % number of discretized fields
x = linspace(0,L,n+1);

% Initial conditions
u0 = .1*ones(1,n-1);
u0 = u0';

dudt = @(t,u) KSd_dirichlet(t,u,v)';

T1 = 100; T2 = 10000; dt = 0.1; N = 400;

% Integrate using matlab ode45
[~,u] = ode15s(dudt,0:dt:T1,u0);
u0 = u(end,:);

[~,u] = ode15s(dudt,0:dt:T2,u0);
u = [zeros(T2/dt+1,1),u,zeros(T2/dt+1,1)];

% % draw contour
% [x,t] = meshgrid(0:L/n:L,20*dt:20*dt:T2);
% [C,h] = contourf(x,t,u(20:20:end,:));
% set(h,'LineColor','none'); colormap(jet); colorbar;
% xlabel('$x$','Interpreter','latex');
% ylabel('$t$','Interpreter','latex');

%%%%%%%%%%%% POD %%%%%%%%%%%%
u_mean = mean(u);
figure(3); hold on; plot(x,u_mean);
u_centered = u - repmat(u_mean,size(u,1),1);
[~,S,V] = svd(u_centered(10:10:end,:));

% plot 4 high energy mode-shapes
figure(4);
for i = 1:4
    subplot(2,2,i);
    plot(x,-V(:,i)); xlim([0 16]);
end

% Energy spectrum
s = diag(S);
Ek = s.^2/sum(s.^2);
Ek_cum = cumsum(Ek);
figure; plot(Ek_cum); xlim([0 50]);
hold on; plot(0.9*ones(50,1),'-.');
line([10 10],[0 1]);

c = u_centered*V;
c_nrm = bsxfun(@rdivide,c,s'/sqrt(10000));        % standardized mode coeffs
% figure; scatter3(c(:,1),c(:,2),c(:,3));

figure; plot(var(c_nrm));

csvwrite('mean.csv',u_mean);
csvwrite('mode_shape.csv',V);
csvwrite('mode_eigvalues.csv',s);
csvwrite('mode_coeff_TS.csv',c);
csvwrite('mode_coeff_TS_normalized.csv',c_nrm);

% calculate derivative set 
% dt = 0.2
dcdt = zeros(size(c));
dcdt(2:end-1,:) = (c(3:end,:) - c(1:end-2,:))/2/dt;
dcdt(1,:) = (c(2,:) - c(1,:))/dt;
dcdt(end,:) = (c(end,:) - c(end-1,:))/dt;

% Get the most energetic modes
rdim = 20;

u_r = c(:,1:rdim); dudt_r = dcdt(:,1:rdim);

[idx_sorted,D] = knnsearch(u_r,u_r(1,:),'K',size(u_r,1));
figure; plot(D);

[u_r,dudt_r] = DataShrink(u_r,5,dudt_r);

% save c0GPRP2585_val.mat u_r dudt_r V u_mean
