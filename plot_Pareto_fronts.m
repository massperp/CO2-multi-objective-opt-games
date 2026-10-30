% Plot the Pareto fronts for all computed coalition structures.
%
%   Reads the results written by main_generate_MOO_results. The settings
%   come from problem_config.m, so they match the optimization run
%   automatically. Only feasible, non-dominated points are plotted.

path_main = pwd;
addpath(genpath(fullfile(path_main, 'Set-Partition')));

% Load settings from config file
cfg = problem_config();


CS_to_plot = 2:cfg.num_CS;   % index 1 is the grand coalition (a single value)

symsize = 30;
fsize   = 15;
marker  = {symsize, 'MarkerFaceColor', '#0072BD', 'MarkerEdgeColor', 'k', 'LineWidth', 0.2};
latex   = {'FontWeight', 'normal', 'Interpreter', 'latex', 'FontSize', fsize};

fig1 = figure('Position', [100, 100, 1300, 360]);
tlo1 = tiledlayout(fig1, 1, numel(CS_to_plot), 'Padding', 'compact');

Obj_PF_all = cell(cfg.num_CS, 1);
for ind_CS = CS_to_plot
    nexttile(tlo1);
    coal_indices = cfg.coal_indices{ind_CS};
    CS = cfg.coal_struct_list{ind_CS};
    num_coal = numel(CS);
    folder_name_sub = results_path(cfg, path_main, ind_CS);

    if ~isfolder(folder_name_sub)
        warning('No results for coalition structure %d in %s; skipping.', ind_CS, folder_name_sub);
        continue
    end

    % Collect decision variables and constraint values of the final populations
    switch cfg.sol_meth
        case 'MOO'
            Dec = load(fullfile(folder_name_sub, 'Dec.mat')).Dec;
            Con = load(fullfile(folder_name_sub, 'Con.mat')).Con;
        case 'WSM'
            wts = load(fullfile(folder_name_sub, 'wts.mat')).wts;
            Dec = []; Con = [];
            for ind_wts = 1:size(wts, 2)
                Dec = [Dec; load(fullfile(folder_name_sub, sprintf('Dec_%d.mat', ind_wts))).Dec]; %#ok<AGROW>
                Con = [Con; load(fullfile(folder_name_sub, sprintf('Con_%d.mat', ind_wts))).Con]; %#ok<AGROW>
            end
    end

    % Injected mass per coalition (Mt), feasible points only
    feasible = all(Con <= 0, 2);
    if ~all(feasible)
        fprintf('CS %d: discarding %d of %d infeasible points.\n', ind_CS, nnz(~feasible), numel(feasible));
    end
    W = zeros(nnz(feasible), num_coal);
    for ind_C = 1:num_coal
        W(:, ind_C) = cfg.rate2Mt * sum(Dec(feasible, coal_indices{ind_C}), 2);
    end

    % Pareto front (maximization of every coalition's injected mass)
    Obj_PF = W(is_nondominated(-W), :);
    Obj_PF_all{ind_CS} = Obj_PF;

    % Axis limits: from the minimum admissible injection to 5% above the max
    lo = cellfun(@(idx) cfg.rate2Mt * sum(cfg.xlow(idx)), coal_indices)';
    hi = 1.05 * max(Obj_PF, [], 1);
    coal_label = @(C) ['\{', strjoin(arrayfun(@(a) sprintf('a_%d', a), C, 'UniformOutput', false), ', '), '\}'];
    labels = arrayfun(@(k) sprintf('$w(C_%d)$ (Mt)', k), 1:num_coal, 'UniformOutput', false);
    ttl = ['$', strjoin(arrayfun(@(k) sprintf('C_%d = %s', k, coal_label(CS{k})), ...
        1:num_coal, 'UniformOutput', false), ', '), '$'];

    if num_coal == 2
        scatter(Obj_PF(:,1), Obj_PF(:,2), marker{:});
        axis([lo(1) hi(1) lo(2) hi(2)]);
    elseif num_coal == 3
        % The surface is a linear interpolation between the computed points.
        xlin = linspace(min(Obj_PF(:,1)), max(Obj_PF(:,1)), 40);
        ylin = linspace(min(Obj_PF(:,2)), max(Obj_PF(:,2)), 40);
        [X, Y] = meshgrid(xlin, ylin);
        Z = griddata(Obj_PF(:,1), Obj_PF(:,2), Obj_PF(:,3), X, Y, 'linear');
        surf(X, Y, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.6); hold on
        scatter3(Obj_PF(:,1), Obj_PF(:,2), Obj_PF(:,3), marker{:});
        zlabel(labels{3}, latex{:});
        view(120, 10);
        axis([lo(1) hi(1) lo(2) hi(2) lo(3) hi(3)]);
        ax = gca; ax.XAxis.TickLabelRotation = 0;
    end
    xlabel(labels{1}, latex{:});
    ylabel(labels{2}, latex{:});
    title(ttl, latex{:});
    set(gca, 'FontSize', fsize);
end
