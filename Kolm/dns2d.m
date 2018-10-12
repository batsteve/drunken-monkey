function [u1, u2]=dns2d(u1,u2,par)
global L1 L2 n1 n2 x y

u1_vec = reshape(u1,n1*n2,1);
u2_vec = reshape(u2,n1*n2,1);

kk = -7:7;
[xx, yy] = ndgrid(2*pi*(1:n1)./n1, 2*pi*(1:n2)/n2);


TT = par.t_start:par.t_step:(par.t_end-par.t_step);
e_diss_mu = zeros(length(TT), 1);
e_inflow_mu = zeros(length(TT), 1);
f_coeff = zeros([length(TT), length(kk), length(kk)]);
% strain tensor eigenvalues?
% McIlhaney beta? (or alpha, or gamma)
% local hessian quantities?
% FTLE?

% time derivative
% moving average

[kF1, kF2] = force(par);
F1 = ifft2(kF1);
F2 = ifft2(kF2);

tic
fprintf('Beginning integration:  initial transients end at t=%0.2f.\n', ....
         par.t_start);

integrand = @(t, y) oderhs(t, y, par);
[~, F]=ode45(integrand,...
                [0, par.t_start],...
                [u1_vec;u2_vec], par.ode45_options );
u1 = reshape(F(end,1:n1*n2),         n2, n1);
u2 = reshape(F(end,n1*n2+1:2*n1*n2), n2, n1);

t = 1;
calc_derived_quantities();

u1_vec = F(end,1:n1*n2)';
u2_vec = F(end,n1*n2+1:2*n1*n2)';


fprintf('Integration reached the end of the initial transient zone.\n');

for t = 2:length(TT)
    fprintf('Reached step %d (time %0.2f) after %0.2f seconds.\n', ...
        t, TT(t), toc);
    
    [~, F]=ode45(integrand,...
                [TT(t-1) TT(t-1)+par.t_step/2 TT(t)],...
                [u1_vec;u2_vec], par.ode45_options );
    u1 = reshape(F(3,1:n1*n2),         n2, n1);
    u2 = reshape(F(3,n1*n2+1:2*n1*n2), n2, n1);

    % save the vel field
    if par.save_vector_fields
        j=round(TT(t)/par.t_step)+1;
        wfile=[par.data_filename_base, sprintf('%4.4i',j)];
        save(wfile, 'u1','u2');
    end
    
    calc_derived_quantities();
    if par.plot_intermediate_vector_fields
        plot_derived_quantities
    end
    
    u1_vec = F(3,1:n1*n2)';
    u2_vec = F(3,n1*n2+1:2*n1*n2)';
end

fprintf('Integration over after %0.2f seconds!\n', toc);

u1 = reshape(F(3,1:n1*n2),         n2, n1);
u2 = reshape(F(3,n1*n2+1:2*n1*n2), n2, n1);

%plot_time_series();

save(par.indicator_filename, 'par', 'TT', 'e_diss_mu', 'f_coeff',...
     'e_inflow_mu', 'u1', 'u2');





function [] = calc_derived_quantities()
%     local_e_kinetic = (u1.^2 + u2.^2);
%     e_kinetic_mu(t) = mean(local_e_kinetic(:));
%     e_kinetic_var(t) = var(local_e_kinetic(:));
%     
    dudx = -1/(2*L1/n1)*([u1(end, :); u1(1:end-1,:)] - [u1(2:end, :); u1(1, :)]);
    dudy = -1/(2*L2/n2)*([u1(:, end), u1(:, 1:end-1)] -[u1(:, 2:end), u1(:, 1)]);
    dvdx = -1/(2*L1/n1)*([u2(end, :); u2(1:end-1,:)] - [u2(2:end, :); u2(1, :)]);
    dvdy = -1/(2*L2/n2)*([u2(:, end), u2(:, 1:end-1)] -[u2(:, 2:end), u2(:, 1)]);
   
    local_e_diss = par.nu*(dudx.^2 + dudy.^2 + dvdx.^2 + dvdy.^2);
    e_diss_mu(t) = mean(local_e_diss(:));

    local_e_influx = (u1.*F1 + u2.*F2);
    e_inflow_mu(t) = mean(local_e_influx(:));


    for kn1 = 1:length(kk)
        for kn2 = 1:length(kk)
            k1 = kk(kn1);
            k2 = kk(kn2);
            e_fac = exp(-1i*(k1.*xx(:) + k2*yy(:)));
            %ff = (u1(:)*k1 + u2(:)*k2)./sqrt(k1^2 + k2^2).* e_fac;
            %f_coeff_1(t, kn1, kn2) = sum(ff(:));
            
            %ff = (u1(:)*k2 - u2(:)*k1)./sqrt(k1^2 + k2^2).* e_fac;
            %f_coeff_2(t, kn1, kn2) = sum(ff(:));
         
            ff = (u1(:) + 1i*u2(:)).* e_fac;
            f_coeff(t, kn1, kn2) = sum(ff(:));
            %f_coeff_4(t, kn1, kn2) = sum(abs(ff(:).^2));
        end
    end
end

function [] = plot_derived_quantities()
    figure(11);
    clf;
    contourf(x,y,local_e_kinetic);
    axis equal tight; colorbar
    title(sprintf('local kinetic energy at t=%0.2f', TT(t)));
    
    figure(12);
    clf;
    contourf(x,y,local_e_diss);
    axis equal tight; colorbar
    title(sprintf('local energy dissipation at t=%0.2f', TT(t)));
    
    figure(13);
    clf;
    contourf(x,y,local_e_influx);
    axis equal tight; colorbar
    title(sprintf('local energy inflow at t=%0.2f', TT(t)));
    
    figure(14);
    clf;
    contourf(x,y,local_vort);
    axis equal tight; colorbar
    title(sprintf('local vorticity at t=%0.2f', TT(t)));
    
    figure(15);
    clf;
    contourf(x,y,local_q);
    axis equal tight; colorbar
    title(sprintf('local Okubu-Weiss at t=%0.2f', TT(t)));
    
    drawnow();
end

    function [] = plot_time_series()
        
        figure(23)
        clf;
        hold on
        kn1 = 1+7;
        kn2 = 2+7;
        YY = (e_diss_mu - min(e_diss_mu(:)))./(max(e_diss_mu(:))- min(e_diss_mu(:)));
        plot(TT, YY);
        %YY = (e_inflow_mu - min(e_inflow_mu))./(max(e_inflow_mu) - min(e_inflow_mu));
        YY = abs(f_coeff_1(:, kn1, kn2))./max(abs(f_coeff_1(:, kn1, kn2)));
        plot(TT, YY);
        legend('dissipated energy', sprintf('%d, %d fourier mode', kk(kn1), kk(kn2)));
        xlabel('T')
        %xlim([100, 500]);
        
        figure(24)
        clf
        scatter(e_diss_mu(10:end), abs(f_coeff_1(1:(end-9), kn1, kn2)), 1);
        
    end

end

