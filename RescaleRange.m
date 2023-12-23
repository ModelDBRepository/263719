function [NewValue]=RescaleRange(OldValue,OldMin,OldMax,NewMin,NewMax)
%% CREDITS
% Created by
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

% Rescaling the values within a range

%% CODE
NewValue = (((OldValue - OldMin) .* (NewMax - NewMin)) ./ (OldMax - OldMin)) + NewMin;

end