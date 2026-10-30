
% Update schedule controls
function [schedule] = rate2scheduleVarRates(u, W, injTimeInts, schedule)

numTimes = length(injTimeInts);
numCtrls = max(schedule.step.control);

for indW = 1:numel(schedule.control(1).W)
    for indCtrls = 1:numCtrls
        schedule.control(indCtrls).W(indW).val = u((indW-1)*numCtrls+indCtrls);    
    end
end
