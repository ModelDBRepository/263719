%% CREDITS
% Created by
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

% RUN script for Levodopa-induced toxicity model
% Included:
% Energy deficiency (soma & terminal)
% Levodopa medication
% SP antagonist therapy
% Glutathione therapy

%% CODE
clc;clear;close all;

dur=10000; % Duration of simulation in milliseconds

ntrials=1; % Number of trials
trial_num=1:1:ntrials;
gpuon=1; % gpuon=1--> Code runs on GPU; gpuon=0--> Code runs on CPU only

wstsn=[0.3]; % STN-->SNc weight
scfa=0.00001; % STN-SNc weight scaling factor
Aapopthr=[0.5]; % Apopotic signal threshold at which cell degenerates
nR=1; %0-autoreceptors 1-without autoreceptors

% Thresholds
caerthr=[2.15e-3]; % Calcium threshold in ER (soma) (mM)
camtthr=[0.0215];  % Calcium threshold in MT (soma) (mM)
rosthr=0.0147;     % ROS threshold in terminal (mM)

% Energy deficiency
ED=[0]; % Percentage of soma/terminal in energy deficiency
peren=ED;  % Percentage of soma in energy deficiency
perenT=ED; % Percentage of terminal in energy deficiency

% Therapy initiation
cl=[25];  % Percentage of soma loss after which therapy was initiated
clT=[25]; % Percentage of terminal loss after which therapy was initiated
Son=0; % Son=1--> Therapy initiated based on soma loss
Ton=0; % Ton=1--> Therapy initiated based on terminal loss
% Son/Ton=0 then no therapy is initiated
if Son==0;cl=0;end
if Ton==0;clT=0;end

% Extracellular levodopa concentration
ldopaN=[36e-4]; % Extracellular levodopa concentration (mM)

% Levodopa therapy
LDOPAon=0; % Levodopa therapy LDOPAon=1--> ON; LDOPAon=0--> OFF
ldopa1=[0]; % Levodopa dosage administrated in soma/terminal (mM)
ldopaS_dose=ldopa1; % Levodopa dosage administrated in soma (mM)
ldopaT_dose=ldopa1; % Levodopa dosage administrated in terminal (mM)

% SP antagonist therapy
SPon=0; % SP antagonist therapy SPon=1--> ON; SPon=0--> OFF
spt_dose=[0]; % SP antagonist dosage

% Glutathione therapy
GSon=0; % Glutathione therapy GSon=1--> ON; GSon=0--> OFF
gst_dose=[0];% >2.5 Glutathione concentration (mM)

wsp=5000;
scfa1=deci2str(scfa);
durr=deci2str(dur/1000);

tic
for km=1:numel(trial_num)
    for i=1:numel(peren)
        for ki=1:numel(perenT)
            for j=1:numel(wstsn)
                for k=1:numel(Aapopthr)
                    for l=1:numel(camtthr)
                        for m=1:numel(cl)
                            for mm=1:numel(clT)
                                for ld=1:numel(ldopa1)
                                    for sp=1:numel(spt_dose)
                                        for gs=1:numel(gst_dose)
                                            
                                            
                                            wwstsn=deci2str(wstsn(j));
                                            peren1=deci2str(peren(i));
                                            perenT1=deci2str(perenT(ki));
                                            apopthr=deci2str(Aapopthr(k));
                                            camtthr1=deci2str(camtthr(l));
                                            cl1=deci2str(cl(m));
                                            clT1=deci2str(clT(mm));
                                            rosthr1=deci2str(rosthr);
                                            wsp1=deci2str(wsp);
                                            ldopaS_dose1=deci2str(ldopaS_dose(ld));
                                            ldopaT_dose1=deci2str(ldopaT_dose(ld));
                                            spt_dose1=deci2str(spt_dose(sp));
                                            ldopaN1=deci2str(ldopaN);
                                            gst_dose1=deci2str(gst_dose(gs));
                                            filename=strcat('GST_LIT_LD_g0-1_m0-15_0-15_LDN',num2str(ldopaN1),'_CT',num2str(camtthr1),'_RT',num2str(rosthr1),'_PS',num2str(peren1),'_PT',num2str(perenT1),'_CL',num2str(cl1),'%_CLT',num2str(clT1),'%_LDS',num2str(ldopaS_dose1),'_LDT',num2str(ldopaT_dose1),'_SPT',num2str(spt_dose1),'_GST',num2str(gst_dose1),'_',num2str(trial_num(km)),'_',num2str(durr),'sec');
                                            
                                            %  filename=strcat('SPT_LIT_LD_g0-1_m0-15_0-15_LDN',num2str(ldopaN1),'_CT',num2str(camtthr1),'_RT',num2str(rosthr1),'_PS',num2str(peren1),'_PT',num2str(perenT1),'_CL',num2str(cl1),'%_CLT',num2str(clT1),'%_LDS',num2str(ldopaS_dose1),'_LDT',num2str(ldopaT_dose1),'_SPT',num2str(spt_dose1),'_',num2str(trial_num(km)),'_',num2str(durr),'sec');
                                            %  filename=strcat('LDT_LIT_LD_g0-1_m0-15_0-15_LDN',num2str(ldopaN1),'_CT',num2str(camtthr1),'_RT',num2str(rosthr1),'_PS',num2str(peren1),'_PT',num2str(perenT1),'_CL',num2str(cl1),'%_CLT',num2str(clT1),'%_LDS',num2str(ldopaS_dose1),'_LDT',num2str(ldopaT_dose1),'_',num2str(trial_num(km)),'_',num2str(durr),'sec');
                                            
                                            disp(filename)
                                            [deda,datp,dedaT,dcdaT,datpT,drosT,indsappT,indsapp,dncS,dncT,dDAT,dIsp,dcaiS,dcaiT,dcamtS,simtime,srnd]=MAIN_LIT_model(dur,peren(i),perenT(ki),wstsn(j),scfa,Aapopthr(k),camtthr(l),rosthr,cl(m),clT(mm),wsp,LDOPAon,ldopaS_dose(ld),ldopaT_dose(ld),Son,Ton,ldopaN,SPon,spt_dose(sp),GSon,gst_dose(gs),gpuon,filename);
                                            parsave_terminal_soma(filename,deda,datp,dedaT,dcdaT,datpT,drosT,indsappT,indsapp,dncS,dncT,dDAT,dIsp,dcaiS,dcaiT,dcamtS,simtime);
                                            disp(filename)
                                            
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
toc