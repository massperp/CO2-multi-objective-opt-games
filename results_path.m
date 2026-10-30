function folder = results_path(cfg, path_main, ind_CS)
%RESULTS_PATH Folder for the results of a run.
%
%   results_path(cfg, path_main)         main results folder
%   results_path(cfg, path_main, ind_CS) subfolder for coalition structure ind_CS
%
%   The folder names include the algorithm parameters, so runs with
%   different settings do not overwrite each other.

    folder = fullfile(path_main, sprintf('results_%s_NW_%d', cfg.sol_meth, cfg.num_wells));
    if nargin < 3
        return
    end
    switch cfg.sol_meth
        case 'MOO'
            name = sprintf('CS_ind_%d_N%d_maxFE%d_SO_init_%d', ind_CS, cfg.N, cfg.maxFE, cfg.maxFE_init);
        case 'WSM'
            name = sprintf('CS_ind_%d_N%d_maxFE%d', ind_CS, cfg.N, cfg.maxFE);
        otherwise
            error('Unknown sol_meth ''%s''; use ''MOO'' or ''WSM''.', cfg.sol_meth);
    end
    folder = fullfile(folder, name);
end
