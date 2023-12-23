function [wlatgpe]= weightcal_gpe(DA)
%% CREDITS
% Created by
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

% Computing GPe lateral connection weight matrix

%% CODE
%----------------GPe_lateral connections-----------%
smax = 0.1002;  %Strength of lateral connections in gpe
rs = 1.6; %Radius of lateral connections in gpe
nlatgpe = 15; % number of laterals in gpe

ssmax = (smax.*(1+43.*DA))./(1-(0.1.*DA.*DA));
% ssmax = (smax*(1+400*DA))/(1+(2.9*DA));
% ssmax = smax.*exp(3.91003.*DA); % expontential one (6-3-18)BEST WORKED ONE
% SS 0.1-->8 DA 0-->1
% ssmax = smax.*exp(7.*DA); %27-8-18

wlatgpe = calclatwts(nlatgpe,ssmax, rs);

end