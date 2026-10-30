
% Evaluate the overburden pressure constraint with input from MRST

function [pressure_constr] = evaluate_pressure_constr(u, P_limit, initState, model, schedule) %, scaling)

% Update schedule to new rates
% Injection time steps: schedule.step.val
[schedule] = rate2scheduleVarRates(u, schedule.control(1).W, schedule.step.val, schedule);

% Run the CO2 storage model
[wellSols, states] = ...
   simulateScheduleAD(initState, ...
                      model, ...
                      schedule, ...
                      'NonLinearSolver', NonLinearSolver('useRelaxation', true));
 
 const_all_times = zeros(numel(states),1); 
 const_rel_all_times = zeros(numel(states),1); 
 for j=1:numel(states)
     const_all_times(j,1) = max(states{j}.pressure - P_limit);
     const_rel_all_times(j,1) = max( (states{j}.pressure - P_limit)./initState.pressure);
 end

 pressure_constr = max(const_rel_all_times);
 