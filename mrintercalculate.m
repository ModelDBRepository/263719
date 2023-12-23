function [Rvalue,Ravg]=mrintercalculate(linear_Stn,linear_Gpe,Nneur,Ttime)

%% CREDITS
% Created by
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

% Computing synchrony between two neural populations

%% CODE
% [m,n]=size(linear_Stn);
Rvalue=[];phistn=zeros(Nneur,Ttime);
phigpe=zeros(Nneur,Ttime);
for neur=1:Nneur
    if isempty(linear_Stn)
        continue
    end
    temptime1=linear_Stn((linear_Stn(:,2)==neur));
    temptime2=linear_Gpe((linear_Gpe(:,2)==neur));
%     temp1=linear_Stn(neur,:);
%     temp2=linear_Gpe(neur,:);
%     temptime1=find(temp1==1);
%     temptime2=find(temp2==1);
j=1;k=1;
for i=1:Ttime % initial transients
       if j<numel(temptime1)
        if i>=temptime1(j) && i<=temptime1(j+1)
            phistn(neur,i)=((2*pi*(i-temptime1(j)))/(temptime1(j+1)-temptime1(j)));
            if i==temptime1(j+1)
                j=j+1;
            end
        end
       end 

    if k<numel(temptime2)
        if i>=temptime2(k) && i<=temptime2(k+1)
            phigpe(neur,i)=((2*pi*(i-temptime2(k)))/(temptime2(k+1)-temptime2(k)));
            if i==temptime2(k+1) 
                k=k+1;
            end
        end
    end
    
end

end

phistn1=phistn(:,:);
phigpe1=phigpe(:,:);

a=sqrt(-1);
tempMstn=sum(phistn1)/numel(phistn1);
tempMgpe=sum(phigpe1)/numel(phigpe1);
M=exp(a*((tempMstn+tempMgpe)/2));
M=M';

sumphistn=(sum(exp(a*phistn1))/neur);
sumphigpe=(sum(exp(a*phigpe1))/neur);
sumphistn=sumphistn';
sumphigpe=sumphigpe';
Rvalue=((sumphistn+sumphigpe)/2)./M;
Ravg=sum(abs(Rvalue))/numel(Rvalue);
% figure(fignumber)
% subplot(1,2,1)
% plot(abs(Rvalue));
% xlabel('time(ms)');
% ylabel('Sync parameter');
% title('Synchronization parameter');

end