function [Rvalue,Ravg]=mrcalculate(linear_S,Nneur,Ntime)

%% CREDITS
% Created by
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

% Computing synchrony

%% CODE
Rvalue=[];phi=[];

phi=3000*ones(Nneur,Ntime-1);
for neur=1:Nneur
        if isempty(linear_S)
            continue
        end
    temptime=linear_S((linear_S(:,2)==neur));
    % temptime =[4    12    21    30    60    78   100   117   126   163   503   652   797   857   940   943];
    j=1;k=1;
    
    while j<numel(temptime)
        %     k=1;
        if temptime(j)==0
            continue
        end
        for i=temptime(j):1:temptime(j+1)-1
            
            phi(neur,i)=(2*pi*(i-temptime(j)))/(temptime(j+1)-temptime(j));
            %                 k=k+1;
        end
        j=j+1;
    end
    
end
phi;
neur;
a=sqrt(-1);
tempM=sum(phi)/numel(phi);
M=exp(a*tempM);
Rvalue=((sum(exp(a*phi))/neur))./M;
Ravg=sum(abs(Rvalue))/numel(Rvalue);


end