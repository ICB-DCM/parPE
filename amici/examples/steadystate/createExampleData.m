function [ output_args ] = createExampleData(  )
    % Create data for steadystate-model parameter estimation example
    
    cd(fileparts(mfilename('fullpath')));
    hdfFile = 'data.h5';
    
    % add example directory to path
    exdir = fullfile(fileparts(which('amiwrap.m')), 'examples', 'example_steadystate');
    addpath(exdir);
    
    % compile the model
    amiwrap('model_steadystate', 'model_steadystate_syms', exdir);
    
    %% model & data dimensions
    model = model_steadystate_syms;
    t = linspace(0,100,20); % timepoints
    t = [t, inf];
    p = [1;0.5;0.4;2;0.1]; % optimization parameters
    assert(numel(p) == numel(model.sym.p));
    % k_original = [0.1,0.4,0.7,1]; % fixed parameters
    
    options = amioption('sensi', 0, ...
                        'maxsteps', 1e4, ...
                        'pscale', 'log10');

    numConditions = 12;
    numK = numel(model.sym.k); % number of constant parameters
    numP = numel(model.sym.p); % number of variable parameters
    numY = numel(model.sym.y); % number of observables
    numT = numel(t); % number of timepoints
    sigmaY = 0.1; % measurement standard deviation
    
    rng(0);
    
    kAll = ones(numConditions, numK);
    
    %% write data
    delete(hdfFile);
    try
        h5create(hdfFile, '/data/k', size(kAll));
    end
    try
        h5create(hdfFile, '/data/ytrue', [numT, numY, numConditions]);
    end
    try
        h5create(hdfFile, '/data/ymeasured', [numT, numY, numConditions]);
    end
    try
        h5create(hdfFile, '/data/sigmay', [numT, numY, numConditions]);
    end
    
    h5write(hdfFile, '/data/k', kAll);
    h5writeatt(hdfFile, '/data/', 'ptrue', log10(p));
    h5writeatt(hdfFile, '/data/', 't', t);
    h5write(hdfFile, '/data/sigmay', sigmaY * ones(numT, numY, numConditions));

    for i = 1:numConditions
        k = kAll(i, :);
        sol = simulate_model_steadystate(t, p, k, [], options);
        h5write(hdfFile, '/data/ytrue', sol.y, [1, 1, i], [numT, numY, 1]);
        h5write(hdfFile, '/data/ymeasured', sol.y * (1 + sigmaY * randn()), [1, 1, i], [numT, numY, 1]);
    end
    
    %% set bounds
    h5create(hdfFile, '/parameters/lowerBound', [numP])
    h5write(hdfFile, '/parameters/lowerBound', -2 * ones(numP, 1))
    h5create(hdfFile, '/parameters/upperBound', [numel(p)])
    h5write(hdfFile, '/parameters/upperBound', 2 * ones(numP, 1))
    
    %% set optimization options
    % create group
    fid = H5F.open(hdfFile, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
    gid = H5G.create(fid, '/optimizationOptions', 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
    H5G.close(gid);
    H5F.close(fid);

    h5writeatt(hdfFile, '/optimizationOptions/', 'optimizer', 0);
    h5writeatt(hdfFile, '/optimizationOptions/', 'maxIter', 20);

end

