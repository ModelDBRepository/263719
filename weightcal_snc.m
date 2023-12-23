function [wlatsnc]= weightcal_snc(DA,rs,nlatsnc)

%% CREDITS
% Created by
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

% Computing SNc lateral connection weight matrix

%% CODE
%----------------SNc_lateral connections-----------%
smax = 0.000001;  %Strength of lateral connections in snc
% rs = 1.6; %Radius of lateral connections in snc
% nlatsnc = 15; % number of laterals in snc

if DA>1.1
    DA=1;
end

if DA<0
    DA=0;
end

ssmax = (smax.*(exp(DA.*4.6055)));

wlatsnc = calclatwts(nlatsnc,ssmax, rs);

end