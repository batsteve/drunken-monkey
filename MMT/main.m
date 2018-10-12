n_sims = 1;

datapath = '../Data/MMT';

if ~exist(datapath, 'dir')
    mkdir(datapath);
end

sim_datapath = sprintf('%s/sim_data', datapath);
if (~exist(sim_datapath, 'dir'))
    mkdir(sim_datapath)
end


for k = 1:n_sims
    fprintf('Starting round %d.\n', k);
    m_par = MMT_Parameters();

    % Integrate MMT NL Equation
    
    filename = sprintf('%s/sim_%d.dat', sim_datapath, k);
    if (~exist(filename, 'file'))
        YY = simulate_MMT_NL_data(m_par);
        save(filename, 'm_par', 'YY');
    end

end