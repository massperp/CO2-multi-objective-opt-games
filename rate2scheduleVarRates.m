
% Written by PP based on example from O. Møyner
function [schedule] = rate2scheduleVarRates(u, W, injTimeInts, schedule)

%%%schedule = simpleSchedule(injTimeInts, 'W', W);

numTimes = length(injTimeInts);
numCtrls = max(schedule.step.control);
%schedule.step.control(1:numTimes) = 1:numTimes;
%schedule.control = repmat(schedule.control, numTimes, 1);
%if nargin < 3
    %%schedule = simpleSchedule(injTimeInts, 'W', W);
    
    % Set a control for every time step
    %%schedule.step.control = 1:numTimes;

    % Repeat the single control struct three times - this duplicates the wells.
    %%schedule.control = repmat(schedule.control, numel(injTimeInts), 1);
%end



% Define three different periods:
% 0.0 -> 0.5 year (control 1)
% 0.5 -> 1.5 year (control 2)
% 1.5 -> 2.0 year (control 3)
% Each timestep has assigned a integer corresponding to a control period.
% % schedule.step.control(time < 0.5*year) = 1;
% % schedule.step.control(time >= 0.5*year & time < 1.5*year) = 2;
% % schedule.step.control(time >= 1.5*year) = 3;

% % disp(schedule.step.control)
% % % Repeat the single control struct three times - this duplicates the wells.
% % schedule.control = repmat(schedule.control, 3, 1);
% % 
% % % Assign different rates to the three injection periods for the injector
% % irate = sum(poreVolume(G, rock))/T;
% % schedule.control(1).W(1).val = 0.25*irate;
% % schedule.control(2).W(1).val = 0.75*irate;
% % schedule.control(3).W(1).val = 0.50*irate;

%%%%%%%%%%%%%


% % for indW = 1:numel(W)
% %     schedule.control(1:numTimes).W(indW).val = u((indW-1)*numTimes+1:indW*numTimes);
% % end
% Need double for loop due to cell indexing, consider transform to vector,
% use vector indices, and transform back



for indW = 1:numel(schedule.control(1).W)

    for indCtrls = 1:numCtrls
        schedule.control(indCtrls).W(indW).val = u((indW-1)*numCtrls+indCtrls);
    
    end
end
