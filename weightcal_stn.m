function [wlatstn]= weightcal_stn(DA)
%% CREDITS
% Created by
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

% Computing STN lateral connection weight matrix

%% CODE
%----------------STN_lateral connections-----------%
smax = 1.3;  %Strength of lateral connections in stn 2.4-dt0.05 1.3-dt0.1
rs = 1.4; %Radius of lateral connections in stn
nlatstn = 11; % number of laterals in stn
% DA=1;
% ssmax = (smax*(1+0.01*DA))/(1+4*DA); % Bursting to single spike syn
% ssmax = (smax.*(1+0.01.*DA))./(1+4.*DA); % Single spike syn to no syn
% ssmax = (smax*(1+20.5*DA))/(1+30*DA);
% ssmax = RescaleRange(DA,0,1,smax,0.01); %5-3-18
ssmax = smax.*exp(-4.87.*DA); % expontential one (27-8-18)BEST WORKED ONE %15-dt0.05 %4.87-dt0.1
% SS 2.4-->0.01 DA 0-->1

wlatstn = calclatwts(nlatstn,ssmax, rs);

end