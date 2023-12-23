function [test2]=sub_sampling_GPU(input1,dt)
%% CREDITS
% Created by
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

% Subsampling of array (GPU enabled)

%% CODE
%     input=dVsnc1;dt=0.005;
input=gpuArray(input1);
[mm,nn]=size(input);
Ntime=(nn)*dt;

test2=zeros(mm,Ntime);%,'gpuArray');
% for tr=1:mm
    
    start=1;stop=1/dt;
    for te=1:Ntime
        
        test2(:,te)=gather(mean(input(:,start:stop),2));
        
        start=start+1/dt;stop=stop+1/dt;
        %         disp(te)
    end
    
% end

% test1=gather(test2);

end