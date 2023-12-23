function ff=deci2str(x)
%% CREDITS
% Created by
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

% Convert decimal point to string

%% CODE
% d=mod(x,1);
f=strcat(num2str(x));
ff=strrep(f,'.','-');

end