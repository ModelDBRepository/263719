function [test2]=sub_sampling(input,dt)
%% CREDITS
% Created by
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

% Subsampling of array

%% CODE

%     input=dVsnc1;dt=0.005;
[mm,nn]=size(input);
Ntime=(nn)*dt;

test2=zeros(mm,Ntime);
for tr=1:mm
    
    start=1;stop=1/dt;
    for te=1:Ntime
        
        test2(tr,te)=mean(input(tr,start:stop));
        
        start=start+1/dt;stop=stop+1/dt;
        %         disp(te)
    end
    
end

end