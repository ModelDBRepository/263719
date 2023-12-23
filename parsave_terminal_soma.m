function parsave_terminal_soma(filename,deda,datp,dedaT,dcdaT,datpT,drosT,indsappT,indsapp,dncS,dncT,dDAT,dIsp,dcaiS,dcaiT,dcamtS,simtime)
%% CREDITS
% Created by
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

% Saving the data

%% CODE
save(filename,'deda','datp','dedaT','dcdaT','datpT','drosT','indsappT','indsapp','dncS','dncT','dDAT','dIsp','dcaiS','dcaiT','dcamtS','simtime','filename','-v7.3')
end