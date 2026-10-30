function cfg = problem_config()
% Shared settings for main_generate_MOO_results and plot_Pareto_fronts.
%
%   cfg = problem_config() returns a struct with the problem setup (agents,
%   wells, controls, bounds), the coalition structures, the evolutionary
%   algorithm parameters, and helpers for the results folders. Both the
%   optimization script and the plotting script use this script.
%
%   Units: the design variables u are well rates in m^3/s at the CO2
%   reference density cfg.rhoCref (kg/m^3). One design variable acts for
%   cfg.inj_period years, so the injected mass in Mt is
%   u * cfg.rate2Mt.


    %% Solution method: 'MOO' (NSGA-II) or 'WSM' (weighted sum method, CSO)
    cfg.sol_meth = 'MOO';

    %% Evolutionary algorithm parameters
    cfg.N_SO       = 20;    % population size, single-objective initialization
    cfg.maxFE_init = 800;   % max function evaluations per initial (end) point
    cfg.N          = 40;    % population size, multi-objective optimization
    cfg.maxFE      = 2000;  % max function evaluations, multi-objective optimization
    cfg.seed       = 1;     % random seed, for reproducible runs

    %% Agents, wells and controls
    cfg.num_agents = 3;
    cfg.num_wells  = 9;
    cfg.inj_time   = 40;    % injection time (years)
    cfg.num_ctrls  = 8;     % number of control periods per well
    assert(mod(cfg.inj_time, cfg.num_ctrls) == 0, ...
        'inj_time (%d) must be a multiple of num_ctrls (%d).', cfg.inj_time, cfg.num_ctrls);
    cfg.inj_period = cfg.inj_time / cfg.num_ctrls;   % years per control period
    cfg.num_vars   = cfg.num_wells * cfg.num_ctrls;

    % Wells operated by each agent
    cfg.well_ind_agents = {[1 2 3 4], [5 6 7], [8 9]};

    %% Bounds on the rates (m^3/s at reference density)
    cfg.xup  = repmat(0.30, 1, cfg.num_vars);
    cfg.xlow = repmat(0.01, 1, cfg.num_vars);
    cfg.ind_wells_late_start = [1 4];   % these wells are shut in the first period
    first_ctrl = (cfg.ind_wells_late_start - 1) * cfg.num_ctrls + 1;
    cfg.xlow(first_ctrl) = 0;
    cfg.xup(first_ctrl)  = 0;

    %% Unit conversion
    cfg.rhoCref = 760;                         % kg/m^3, must match make_Bjarmeland_model_inj_only
    cfg.year_s  = 365.2425 * 24 * 3600;        % s, same definition as MRST's year()
    cfg.rate2Mt = cfg.inj_period * cfg.year_s * cfg.rhoCref / 1e9;  % (m^3/s per period) -> Mt

    %% Design-variable indices
    % Variables are ordered well by well: well w owns columns
    % (w-1)*num_ctrls+1 : w*num_ctrls (see rate2scheduleVarRates).
    %ctrl_ind_well = reshape(1:cfg.num_vars, cfg.num_ctrls, cfg.num_wells)';
    %cfg.design_ind_agents = cellfun(@(w) reshape(ctrl_ind_well(w,:)', 1, []), ...
    %    cfg.well_ind_agents, 'UniformOutput', false);


    % Array to match controls with design variables pertaining to particular
    % wells
    ctrl_ind_well = zeros(cfg.num_wells, cfg.num_ctrls);
    for ind_wells = 1:cfg.num_wells
        ctrl_ind_well(ind_wells,:) = cfg.num_ctrls*(ind_wells-1)+1:cfg.num_ctrls*ind_wells;
    end

    cfg.design_ind_agents = cell(cfg.num_agents,1);
    for ind_agents = 1:cfg.num_agents
        for ind_wells = 1:length(cfg.well_ind_agents{ind_agents})
            cfg.design_ind_agents{ind_agents} = [cfg.design_ind_agents{ind_agents} ctrl_ind_well(cfg.well_ind_agents{ind_agents}(ind_wells),:)];
        end
    end

    %% Coalition structures (all set partitions of the agents)
    cfg.coal_struct_list = {};
    for n = 1:cfg.num_agents
        cfg.coal_struct_list = [cfg.coal_struct_list; SetPartition(cfg.num_agents, n)];
    end
    cfg.num_CS = numel(cfg.coal_struct_list);

    % Design-variable indices of each coalition, for each coalition structure:
    % cfg.coal_indices{ind_CS}{ind_C}
    cfg.coal_indices = cell(cfg.num_CS, 1);
    for ind_CS = 1:cfg.num_CS
        CS = cfg.coal_struct_list{ind_CS};
        cfg.coal_indices{ind_CS} = cellfun(@(agents) [cfg.design_ind_agents{agents}]', ...
            CS, 'UniformOutput', false);
    end
end






