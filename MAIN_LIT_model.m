%% Hybrid network model of excitotoxicity
function [deda,datp,dedaT,dcdaT,datpT,drosT,indsappT,indsapp,dncS,dncT,dDAT,dIsp,dcaiS,dcaiT,dcamtS,simtime,srnd]=MAIN_LIT_model(durr,peren,perenT,wstsn,scfa,apopthr,camtthr,rosthr,cl,clT,wsp,LDOPAon,ldopaS_dose,ldopaT_dose,Son,Ton,ldopaN,SPon,spt_dose,GSon,gst_dose,gpuon,filename)
%% CREDITS
% Created by
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

% LITERATURE
% SNc with ATP dynamics (Francis et.al., 2013)
% Dopamine synthesis, storage, release, metabolism and terminal autoreceptors (Bravo, 2012)
% Ca2+ induced apoptosis (Hong et.al., 2012)
% Calcium-induced calcium release (Marhl et.al., 2000)
% Energy Metabolism (Cloutier & Wellstead, 2010)
% STN, GPe - (Alekhya et.al., 2015)

%%
%Created on Wed Oct 18 22:00:48 2018

%@author: Vignan (IIT-Madras)
%%
tic

sLD=ldopaN;
sLDT=ldopaN;

forstart=1;

wctr_msn1=100;
wctr_msn2=100;
w_msn1_msn2=500; %500
w_msn1_snc=0.5;% = 1
w_msn2_snc=0.5;% = 1
% wsp=5000;  %1-SP connectivity 0- no SP connectvity = 1
noDAT=1;
DA_SP=01;
pdon=0;
pd=1;%STN laterals
scfamsn=4.1457e-06;

% Time parameters & Random seeding
dt=0.1;
taustn=dt;
taugpe=dt;
tspan=dt:dt:durr;
Ttime=numel(tspan);
srnd = rng;

%%
%STN
nSTN=32; % (nSTNxnSTN network size)
Mstn=nSTN;
Nstn=nSTN;
Pstn=Mstn*Nstn;

%SNc
nSNc=8; % (nSNcxnSNc network size)
Msnc=nSNc;
Nsnc=nSNc;
Psnc=Msnc*Nsnc;

%GPe
nGPe=32; % (nGPexnGPe network size)
Mgpe=nGPe;
Ngpe=nGPe;
Pgpe=Mgpe*Ngpe;

% Neuron properties
%STN
astn=0.005; % (1/ms)
bstn=0.265; % (1/mV)
cstn=-65; % (mV)
dstn=1.5; % 2-Thibeault2013
vpeak_stn=30;

%GPe
agpe=0.1; % (1/ms)
bgpe=0.2; % (1/mV)
cgpe=-65; % (mV)
dgpe=2;
vpeak_gpe=30;

% Membrane capacitances
Cstn = 1; %(microF)
Cgpe = 1; %(microF)

% V, U initialization
%STN
Vstn = -62.5.*(rand(Mstn,Nstn)-0.5.*ones(Mstn,Nstn));
Ustn = ((-15)-(-5)).*rand(Mstn,Nstn) + (-5);

%GPe
Vgpe = -53.67.*(rand(Mgpe,Ngpe)-0.5.*ones(Mgpe,Ngpe));
Ugpe = ((-15)-(-5)).*rand(Mgpe,Ngpe) + (-5);

% Currents initilization
%STN
Istn=3*ones(size(Vstn)); % 10 Hz->1.9
stn_zeros = zeros(size(Vstn));
stncurr_spk = stn_zeros;

%GPe
Igpe=4.25*ones(size(Vgpe));
gpe_zeros = zeros(size(Vgpe));
gpecurr_spk = gpe_zeros;

% psp variable initilization
%STN
h_nmdastn=stn_zeros;
h_ampastn=stn_zeros;
h_gs = stn_zeros;

%GPe
h_nmdagpe=gpe_zeros;
h_ampagpe=gpe_zeros;
hlat_gaba_gpe=gpe_zeros;

% decay constants(ms)
taunmda=160;
tauampa=6;
taugaba=4;

% dt/T in PSP
lam_nmda = dt/taunmda;
lam_ampa = dt/tauampa;
lam_gaba= dt/taugaba;

mg0=1; % magnesium conc.

% RMP of receptors
Enmda = 0;
Eampa = 0;
Egaba = -60;

xmin=-55;xmax=-45;
V_sncinit = xmin+rand(Msnc,Nsnc)*(xmax-xmin);
Ca_iinit = 0.000188.*ones(Msnc,Nsnc);
Na_iinit = 4.6876.*ones(Msnc,Nsnc);
K_iinit = 126.05893.*ones(Msnc,Nsnc);
Calbinit = 0.0026.*ones(Msnc,Nsnc);
Caminit = 0.0222.*ones(Msnc,Nsnc);
m_calinit = 0.006271.*ones(Msnc,Nsnc);
m_nainit = 0.0952.*ones(Msnc,Nsnc);
h_nainit = 0.1848.*ones(Msnc,Nsnc);
O_hcninit = 0.003.*ones(Msnc,Nsnc);
m_kdrinit = 0.0932.*ones(Msnc,Nsnc);
y_pcinit = 0.483.*ones(Msnc,Nsnc);
y_nkinit = 0.6213.*ones(Msnc,Nsnc);
Ca_erinit = 1.0*0.001.*ones(Msnc,Nsnc); %mM
Ca_mtinit = 0.4*0.001.*ones(Msnc,Nsnc); % mM % 0.4e-3
cdainit = 1e-4.*ones(Msnc,Nsnc); %mM%1e-4
vdainit = 500.*ones(Msnc,Nsnc); %mM 500
edainit = 4e-6.*ones(Msnc,Nsnc); %mM
Iextinit=0.*ones(Msnc,Nsnc); %Iext
ATPusedinit=0.*ones(Msnc,Nsnc);
calinit=1.*ones(Msnc,Nsnc); %mM
cai_calinit=0.*ones(Msnc,Nsnc); %mM
cal_actinit=0.*ones(Msnc,Nsnc); %mM
casp12init=1.*ones(Msnc,Nsnc); %mM
cal_act_casp12init=0.*ones(Msnc,Nsnc); %mM
casp12_actinit=0.*ones(Msnc,Nsnc); %mM
casp9init=1.*ones(Msnc,Nsnc); %mM
casp12_act_casp9init=0.*ones(Msnc,Nsnc); %mM
casp9_actinit=0.*ones(Msnc,Nsnc); %mM
casp3init=1.*ones(Msnc,Nsnc); %mM
casp9_act_casp3init=0.*ones(Msnc,Nsnc); %mM
casp3_actinit=0.*ones(Msnc,Nsnc); %mM
xmin1=0;xmax1=0.1;
apopinit = xmin1+rand(Msnc,Nsnc)*(xmax1-xmin1);
ROS_mitinit=0.*ones(Msnc,Nsnc);
PTP_mit_actinit=0.*ones(Msnc,Nsnc);
Cytc_mitinit=1.*ones(Msnc,Nsnc);
Cytcinit=0.*ones(Msnc,Nsnc);
Cytc_casp9init=0.*ones(Msnc,Nsnc);
IAPinit=0.*ones(Msnc,Nsnc);
casp9_act_IAPinit=0.*ones(Msnc,Nsnc);
casp3_act_IAPinit=0.*ones(Msnc,Nsnc);
F6Pinit = 0.175883476634895.*ones(Msnc,Nsnc);%0.2
F26Pinit = 0.002191750879602.*ones(Msnc,Nsnc);%0.001
GAPinit = 0.082507126186107.*ones(Msnc,Nsnc);%0.0405
PYRinit = 0.123910489378719.*ones(Msnc,Nsnc);%0.1
LACinit = 0.598605032933119.*ones(Msnc,Nsnc);%0.5
ATPinit = 2.395615876085214.*ones(Msnc,Nsnc);%2.402
PCrinit = 18.044071098085976.*ones(Msnc,Nsnc);%18.14
ROSinit=0.0010.*ones(Msnc,Nsnc);%mM
ASYNinit=0.1000.*ones(Msnc,Nsnc);%mM
ASYNAinit=0.0010.*ones(Msnc,Nsnc);%mM
ASYNTinit=9.9997e-06.*ones(Msnc,Nsnc);%mM
ASYNGinit=1.1714e-13.*ones(Msnc,Nsnc);%mM
LBinit=5.0403e-84.*ones(Msnc,Nsnc);%mM
LDOPAinit=3.6000e-04.*ones(Msnc,Nsnc);%mM
NADPHinit=0.2500.*ones(Msnc,Nsnc);%mM
GSHinit=2.5000.*ones(Msnc,Nsnc);%mM
Ssnc=zeros(Msnc,Nsnc);

ROS_mit=ROS_mitinit;
PTP_mit_act=PTP_mit_actinit;
Cytc_mit=Cytc_mitinit;
Cytc=Cytcinit;
Cytc_casp9=Cytc_casp9init;
IAP=IAPinit;
casp9_act_IAP=casp9_act_IAPinit;
casp3_act_IAP=casp3_act_IAPinit;
F6P=F6Pinit;
F26P=F26Pinit;
GAP=GAPinit;
PYR=PYRinit;
LAC=LACinit;
ATP=ATPinit;
PCr=PCrinit;

V_snc=V_sncinit;m_cal=m_calinit;m_na=m_nainit;
h_na=h_nainit;O_hcn=O_hcninit;Calb=Calbinit;
Cam=Caminit;y_nk=y_nkinit;y_pc=y_pcinit;m_kdr=m_kdrinit;
K_i=K_iinit;Na_i=Na_iinit;Ca_i=Ca_iinit;Ca_er=Ca_erinit;Ca_mt=Ca_mtinit;

cda=cdainit;vda=vdainit;
eda=edainit;ATPused=ATPusedinit;
cal=calinit;cai_cal=cai_calinit;cal_act=cal_actinit;casp12=casp12init;
cal_act_casp12=cal_act_casp12init;casp12_act=casp12_actinit;casp9=casp9init;
casp12_act_casp9=casp12_act_casp9init;casp9_act=casp9_actinit;casp3=casp3init;
casp9_act_casp3=casp9_act_casp3init;casp3_act=casp3_actinit;apop=apopinit;

ROS=ROSinit;
ASYN=ASYNinit;
ASYNA=ASYNAinit;
ASYNT=ASYNTinit;
ASYNG=ASYNGinit;
LB=LBinit;
LDOPA=LDOPAinit;
NADPH=NADPHinit;
GSH=GSHinit;

sim_mM=1e-3;
sim_mM_msec=1e-3/3.6e6;
sim_msec=1/3.6e6;
sim_msec_mM=1/((1e-3)*(3.6e6));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
%Soma
R = 8314.472; %mJ/mol. K
T = 310.15; %K
F = 96485.30929; %coulomb/mol.
Ca_o = 1.8; %mM
Na_o = 137; %mM
K_o = 5.4; %mM
vol_pmu = 5; %pl
fr_cyt = 0.5;
C_sp = 0.9e6; %pF/cm2
SVR_pmu = 1.6667e4; %1/cm
% ATP = 0.411; %mM 0.411 - bursting
Calbtot = 0.005; % mM
Camtot = 0.0235; % mM
kcal_1 = 10; %1/mM.ms
kcal_2 = 2e-3; %1/ms
kcam_cd = 0.003; %1/ms
kcam_nd = 3; %1/ms
g_cal = 2101.2; %pA/mM
g_na = 907.68; %pA/mM
A_mna = 1.9651; %1/ms
B_mna = 0.0424; %1/ms
A_hna = 9.566e-5; %1/ms
B_hna = 0.5296; %1/ms
za_mna = 1.7127;
zb_mna = 1.5581;
za_hna = -2.4317;
zb_hna = -1.1868;
g_nalk = 0.0053; %pA/mM
g_nahcn = 51.1; %pA/mM
cAMP = 1e-5; %mM
g_ksk = 2.2515; %pA/mM
g_kdr = 31.237; %nS
g_kir = 13.816; %nS
k_2pc = 0.001; %1/ms
k_3pc = 0.001; %1/ms
k_4pc = 1; %1/ms
K_pco = 2; %mM
k_pmca = 2.233;
dell = 0.35;
k_xm = 0.0166; %pA
k_2nk = 0.04; %1/ms
k_3nk = 0.01; %1/ms
k_4nk = 0.165; %1/ms
K_nknai = 4.05; %mM
K_nknao = 69.8; %mM
K_nkki = 32.88; %mM
K_nkko = 0.258; %mM
k_nk = 1085.7; %pA
V_tau = (R*T)./F;
vol_cyt = fr_cyt*vol_pmu;
P_c = 1.00000./(1.00000+cAMP./0.00116300);
P_o = 1.00000./(1.00000+cAMP./1.45000e-05);
P_E2Spc = 1.00000./(1.00000+K_pco./Ca_o);
A_pmu = (SVR_pmu.*vol_pmu.*0.00100000.*0.00100000.*0.00100000)./1.00000;
P_E2pc = 1.00000-P_E2Spc;
beta_pc = k_2pc.*P_E2Spc+k_4pc.*P_E2pc;

%CICR
% ER
rho_er = 0.01; % rho_er
beta_er = 0.0025; %beta_er
k_pump = 20/1000; % 1/ms
k_ch = 3000/1000; %1/ms
K1 = 5*0.001; %mM
k_leak = 0.05/1000; %1/ms

% Mito
rho_mt = 0.01; % rho_mt
beta_mt = 0.0025; % beta_mt
k_in = 300*0.001/1000; % mM/ms
K2 = 0.8*0.001; % mM
k_out = 125/1000; % 1/ms
k_m = 0.00625/1000; % 1/ms
K3 = 5*0.001; % mM

% DA terminal (Tello-Bravo (2012))
krel = 0.031; %(mM) 0.055
psi = 17.4391793; %(mM/ms)
% nRRP = 10; % :RANGE ~ 10-30
Veda_max = 1e-6; %(mM/ms)
Keda = 3e-5; %(mM)
kcomt = 0.0083511; %(1/ms)
%vda = 500; %(mM)
vdao = 500; %(mM)
vdas = 1e-2; %(mM)
dara = 5e-5; %(mM)
dars = 1e-2; %(mM)
Vsynt_max = 250e-5;%(mM/ms)%30.2e-6 %25e-6
Ksynt = 35e-4; %(mM)
Ktyr = 46e-3; %(mM)
TYR = 126e-3; %(mM)
Kicda = 11e-2; %(mM)
Kieda = 46e-3; %(mM)
Vcda_max = 0.035*133.33e-6; %(mM/ms)%133.33e-6*0.03
Kcda = 238e-4; %(mM)
kmao = 0.00016; %(1/ms)

% Apoptosis pathway (Hong et.al., (2012))
k3f=1; % (muM*sec)-1
k3b=1/1e3; % (sec)-1
k4f=1/1e3; % (sec)-1
k5f=1; % (muM*sec)-1
k5b=1/1e3; % (sec)-1
k6f=1/1e3; % (sec)-1
k7f=10; % (muM*sec)-1
k7b=0.5/1e3; % (sec)-1
k8f=1/1e3; % (sec)-1
k9f=10; % (muM*sec)-1
k9b=0.5/1e3; % (sec)-1
k10f=0.1/1e3; % (sec)-1
k11f=1; % (muM*sec)-1

k29f=0.5; % (mM*msec)-1
k30f=0.5; % (mM*msec)-1
k31f=1; % (mM*msec)-1
k27f=1; % (mM*msec)-1
k27b=1/1e3; % (msec)-1
k28f=1/1e3; % (msec)-1
k12f=5; % (mM*msec)-1
k12b=0.0035/1e3; % (msec)-1
k13f=5; % (mM*msec)-1
k13b=0.0035/1e3; % (msec)-1

Mit=1;
Sig_ers=0;%0.0001;
Sig_mts=0;%0.0001;
PTP_mit=1;

% Energy Metabolism
GLCe=1;%mM
Vmax_hk = 2.5/1000;%mM/ms
Km_ATP_hk = 0.5;%mM
KI_F6P = 0.068;%mM
Vmax_pfk = 3.85/1000;%mM/ms
Km_ATP_pfk = 0.05;%mM
Km_F6P_pfk = 0.18;%mM
Km_F26P_pfk = 0.01;%mM
Vmaxf_pfk2 = 2e-04/1000;%mM/ms
Vmaxr_pfk2 = 1.036e-04/1000;%mM/ms
Km_ATP_pfk2 = 0.05;%mM
Km_F6P_pfk2 = 0.01;%mM
Km_F26P_pfk2 = 0.0001;%mM
Vmax_pk = 5/1000;%mM/ms
Km_ADP_pk = 0.005;%mM
Km_GAP_pk = 0.4;%mM
Vmax_op = 1/1000;%mM/ms
Km_ADP_op = 0.005;%mM
Km_PYR_op = 0.5;%mM
kf_ldh = 12.5/1000;%1/ms
kr_ldh = 2.5355/1000;%1/ms
kf_ck = 3/1000;%1/mM.ms
kr_ck = 1.26/1000;%1/mM.ms
PCr_tot = 20;%mM
Vmax_ATPase = 0.9355/1000;%mM/ms
Km_ATP = 0.5;%mM
Vlac_0 = 0.355/1000;%mM/ms
K_lac_eff = 0.71/1000;%1/ms
K_lac = 0.641;
ANP = 2.51;%mM
Q_adk = 0.92;
nATP = 0.4;
KI_ATP = 1;%mM
nAMP = 0.5;
Ka_AMP = 0.05;%mM
Kamp_pfk2 = 0.005;%mM
nh_amp = 2;
beta_ldh_ros=0.25;
Kldh_ros=10*sim_mM;%muM

kf_gr=0.65.*sim_msec_mM;%1./mM.ms
kr_gr=1.25e-3.*sim_msec_mM;%1./mM.ms
GSH_tot=2500.*sim_mM;%mM
GSH_totT=2500.*sim_mM;%mM
NADPH_tot=250.*sim_mM;%mM
Vmax_ppp=1.43e6.*sim_mM_msec;%mM./ms
Ki_nadph=20;

% PD pathology pathways
eta_op_max=0.995;
beta_op_asyn=0.08;
Kasyn=8.5.*sim_mM; %mM
Kros_cat=235.*1.*sim_msec; %1./ms
Vros_ex=0;%0.1
Kros_dopa=1500.*sim_msec_mM;%1./mM.ms
Kros_dox=0.27.*sim_msec;%1./ms
Kasyn_syn=50.*sim_mM_msec;%mM./ms
Kasyn_ox=7e-5.*sim_msec_mM;%1./mM.ms
Kasyn_to=0.5.*sim_msec;%1./ms
Krasyn_agg=7.5e-4.*sim_msec;%1./ms
Kasyn_agg=7.5.*sim_mM;%mM
Ub_tot=10.5.*sim_mM;%mM
Kasyn_tag=2.75e-7.*sim_msec_mM;%1./mM.ms
Krasyn_prt=7.5e-4.*sim_msec;%1./ms
Kasyn_prt=5.*sim_mM;%mM
beta_asyn_prt=0.25;
Kasyn_lyso=7.5e-5.*sim_msec;%1./ms
Krasyn_lb=7.5e-5.*sim_msec;%1./ms
Kasyn_lb=5.*sim_mM;%
Kros_da=1500*sim_msec_mM;%1/mM.ms
Kda=8.5*sim_mM; %mM

% LDOPA uptake
%Reed et al(2012)
Vaadc_max = 35*2.78e-6;% (mM./ms)
Kaadc = 0.13;% (mM)
Vtran_max = 5.11e-7;% (mM./ms)5.94e-8
% sLD = 36e-3;% (mM)
% sLD = 0;% (mM)
sTYR = 63e-3;% (mM)
sTRP = 82e-3;% (mM)
Ksld = 32e-3;% (mM)
Kstyr = 64e-3;% (mM)
Kstrp = 15e-3;% (mM)
% sLD=3.63685e-3;%35.39059803e-3; %LDOPA=36e-5 Org-36e-3 %mM
vAADC_Km = 0.13;
vAADC_Vmax = 2.78e-6;
eta_nrrp_max=0.995;
beta_nrrp_asyn=0.08;

%% SNc Terminal

%SNc Terminal
nSNcT=32; % (nSNcTxnSNcT network size)
MsncT=nSNcT;
NsncT=nSNcT;
PsncT=MsncT*NsncT;

cdainitT = 1e-4.*ones(MsncT,NsncT); %mM%1e-4
vdainitT = 500.*ones(MsncT,NsncT); %mM 500
edainitT = 4e-6.*ones(MsncT,NsncT); %mM
F6PinitT = 0.175883476634895.*ones(MsncT,NsncT);%0.2
F26PinitT = 0.002191750879602.*ones(MsncT,NsncT);%0.001
GAPinitT = 0.082507126186107.*ones(MsncT,NsncT);%0.0405
PYRinitT = 0.123910489378719.*ones(MsncT,NsncT);%0.1
LACinitT = 0.598605032933119.*ones(MsncT,NsncT);%0.5
ATPinitT = 2.395615876085214.*ones(MsncT,NsncT);%2.402
PCrinitT = 18.044071098085976.*ones(MsncT,NsncT);%18.14
ROSinitT=0.0010.*ones(MsncT,NsncT);%mM
ASYNinitT=0.1000.*ones(MsncT,NsncT);%mM
ASYNAinitT=0.0010.*ones(MsncT,NsncT);%mM
ASYNTinitT=9.9997e-06.*ones(MsncT,NsncT);%mM
ASYNGinitT=1.1714e-13.*ones(MsncT,NsncT);%mM
LBinitT=5.0403e-84.*ones(MsncT,NsncT);%mM
LDOPAinitT=3.6000e-04.*ones(MsncT,NsncT);%mM
NADPHinitT=0.2500.*ones(MsncT,NsncT);%mM
GSHinitT=2.5000.*ones(MsncT,NsncT);%mM

F6PT=F6PinitT;
F26PT=F26PinitT;
GAPT=GAPinitT;
PYRT=PYRinitT;
LACT=LACinitT;
ATPT=ATPinitT;
PCrT=PCrinitT;
cdaT=cdainitT;
vdaT=vdainitT;
edaT=edainitT;
ROST=ROSinitT;
ASYN_T=ASYNinitT;
ASYNAT=ASYNAinitT;
ASYNTT=ASYNTinitT;
ASYNGT=ASYNGinitT;
LBT=LBinitT;
LDOPAT=LDOPAinitT;
NADPHT=NADPHinitT;
GSHT=GSHinitT;
% sLDT=3.63685e-3;%35.39059803e-3; %LDOPA=36e-5 Org-36e-3 %mM
Vros_exT=0;%0.1
adaT=1;

%SNc soma to terminal projections (no. of Soma (1) to no. (proj.*proj) terminal)
proj=4;
nproj=proj*proj;
% idx_soma = randsample(1:Psnc,Psnc);
idx_terminal = randsample(1:PsncT,PsncT);

%% MSN

%msn
nmsn=32; % (nmsnxnmsn network size)
Mmsn=nmsn;
Nmsn=nmsn;
Pmsn=Mmsn*Nmsn;

%%%%%%%%%%%%%%%%%%%%%%%%%%   msn1 Parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
amsn1=0.01	; % (1/ms)
bmsn1=-20; % (1/mV)
cmsn1=-55; % (mV)
dmsn1=91;
vpeak_msn1=40;%(mV)
kmsn1=1;
Kmsn1=0.0289;
Lmsn1=0.331;
Vmsn1r=-80;%(mV)
Vmsn1t=-29.7;%(mV)
% alpha_msn1=0.032;
% phi1=0;
% phi2=0;
% taumsn1 = 0.1

%%%%%%%%%%%%%%%%%%%%%%%%%%   msn2 Parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
amsn2=0.01	; % (1/ms) % Cullen2015
bmsn2=-20; % (1/mV)
cmsn2=-55; % (mV)
dmsn2=91;
vpeak_msn2=40;%(mV)
kmsn2=1;
Kmsn2=0.0289;
Lmsn2=0.331;
Vmsn2r=-80;%(mV)
Vmsn2t=-29.7;%(mV)
% alpha_msn2=0.032;
phi1=0;
% phi2=0;
% taumsn2 = 0.1;

taumsn1=dt;
taumsn2=dt;
tauctr=dt;


%msn1 model
Vmsn1r=Vmsn1r.*(1+(Kmsn1.*phi1));
dmsn1=dmsn1.*(1-(Lmsn1.*phi1));
% kmsn1=kmsn1.*(1-(alpha_msn1.*phi2));

%msn2 model
Vmsn2r=Vmsn2r.*(1+(Kmsn2.*phi1));
dmsn2=dmsn2.*(1-(Lmsn2.*phi1));
% kmsn2=kmsn2.*(1-(alpha_msn2.*phi2));

% Membrane capacitances
Cmsn1 = 15.2; %(pF)
Cmsn2 = 15.2; %(pF)

% V, U initialization
%------MSN1------%
% Vmsn1 = -62.0.*(rand(Mmsn,Nmsn)-0.5.*ones(Mmsn,Nmsn));
% Umsn1=-311.81.*ones(size(Vmsn1));
Vmsn1 = -62.0.*(rand(Mmsn,Nmsn)-0.5.*ones(Mmsn,Nmsn));
Umsn1 = ((-300)-(-200)).*rand(Mmsn,Nmsn) + (-200);
%------MSN2------%
% Vmsn2 = -62.0.*(rand(Mmsn,Nmsn)-0.5.*ones(Mmsn,Nmsn));
% Umsn2=-311.81.*ones(size(Vmsn2));
Vmsn2 = -62.0.*(rand(Mmsn,Nmsn)-0.5.*ones(Mmsn,Nmsn));
Umsn2 = ((-300)-(-200)).*rand(Mmsn,Nmsn) + (-200);
beta1=0.5;

% msn1 current initilization
%msn1
% Iext=270*(1+beta1*phi1);
% Imsn1=Iext*ones(size(Vmsn1));
% Imsn1=gaussDistrMatrix(Mmsn1,Nmsn1,9,0.2); % std = 0.2 from (Lindahl-2013)
% Imsn1=4.5*randn(size(Vmsn1));
msn1_zeros = zeros(size(Vmsn1));
msn1curr_spk = msn1_zeros;
% spkhstmsn1=[];
% msn1_firings=[];

% msn2 current initilization
%msn2
% Iext=270*(1+beta1*phi1);
% Imsn2=Iext*ones(size(Vmsn2));
% Imsn2=gaussDistrMatrix(Mmsn2,Nmsn2,9,0.2); % std = 0.2 from (Lindahl-2013)
% Imsn2=4.5*randn(size(Vmsn2));
msn2_zeros = zeros(size(Vmsn2));
msn2curr_spk = msn2_zeros;
% spkhstmsn2=[];
% msn2_firings=[];

% psp variable initilization
%MSN------%
% hlat_gaba_snc=0;
hlat_gaba_msn=0;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Cortex  %%%%%%%%%%%%%%%%%%%%%%%%
Iextctr=100;%lowest is fires = 52

%CTR
nCTR=32; % (nCTRxnCTR network size)
Mctr=nCTR;
Nctr=nCTR;
Pctr=Mctr*Nctr;

% Neuron properties
%CTR
actr=0.03; % (1/ms)
bctr=-2; % (1/mV)
cctr=-50; % (mV)
dctr=100; % 2-Thibeault2013
vpeak_ctr=35;

% Membrane capacitances
Cctr = 100; %(microF)

% V, U initialization
%CTR
Vctr = -62.5.*(rand(Mctr,Nctr)-0.5.*ones(Mctr,Nctr));
Uctr = ((50)-(-250)).*rand(Mctr,Nctr) + (-250);

% Currents initilization
%CTR
Ictr=Iextctr*ones(size(Vctr)); % 10 Hz->1.9; 3 for 14Hz
ctr_zeros = zeros(size(Vctr));
% ctrcurr_spk = ctr_zeros;
% ctr_firings=[];
% dVctr1=zeros(1,Ttime);

% Ictr_msn1=[];
% Ictr_msn2=[];

Vmsn11r=Vmsn1r;
Vmsn22r=Vmsn2r;
dmsn11=dmsn1;
dmsn22=dmsn2;
%----------------------Array-----------------%
% dVmsn1=[];
% dVsnc1=[];
% Iz_mss=[];
h_nmda_msn1=0;
h_nmda_msn2=0;
h_ampa_msn1=0;
h_ampa_msn2=0;

%%
% snc_firings=[];snc_firings2=[];
DA=0;rss=1.6;nlat=5;

%STN_SNc projections (no. of STN (projXproj) to no. (1) of SNc)
proj=4;
nproj=proj*proj;
idx = randsample(1:Pstn,Pstn);

%SNc ro MSN projections (no. of SNc terminal (projXproj) to no. (1) of MSN)
proj_DA_MSN=20;idx_DA_MSN=zeros(Pmsn,proj_DA_MSN);
for mm=1:Pmsn
    idx_DA_MSN(mm,:)=randperm(PsncT,proj_DA_MSN);
end

%MSN ro SNc projections (no. of MSN (projXproj) to no. (1) of SNc)
proj_MSN_SNc=200;idx_MSN_SNc=zeros(Psnc,proj_MSN_SNc);
for mm=1:Psnc
    idx_MSN_SNc(mm,:)=randperm(Pmsn,proj_MSN_SNc);
end

%SNc
h_nmdasnc=zeros(Msnc,Nsnc);
h_ampasnc=zeros(Msnc,Nsnc);
hlat_msn1_snc=zeros(Msnc,Nsnc);
hlat_msn2_snc=zeros(Msnc,Nsnc);

% STN-SNc connections
% stsn=[1];
wstnsnc_matrix=wstsn.*normrnd(wstsn,0.5,8,8);%sd=0.1

% Effect of DA on post-synaptic currents
cd2=0.1;CD2=0.1;

wsg1=1;wgs1=20;

da=0.5;DAA=1;
wlatgpe = weightcal_gpe(da);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Stimulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Istim=0; % Phasic bursting Istim=0.00001
del=1000/dt;
dur=200/dt;
% sigg1=0;sigg2=0;phier=0;phimt=0;
% enfatra=1;
% enfamit=1;
eda1=0;

% dDA=[];
% dcai=[];datpused=[];dapop=[];
% deda=[];
% dV_snc=[];dcaer=[];dcamt=[];dcam=[];dcalb=[];dros_mit=[];

% dIstnsnc=[];dIsncsnc=[];dIstnstn=[];dIgpestn=[];dIgpegpe=[];dIstngpe=[];dV_snc=[];datp=[];

% sncstart=500;
count=1;
stt=1000;
samp_start=stt/dt;
% yeda1=zeros(Psnc,samp_start);
% dapop1=zeros(Psnc,samp_start);
% dnc1=zeros(1,samp_start);
% dDA1=zeros(1,samp_start);
% dros_mit1=zeros(Psnc,samp_start);
% dcai1=zeros(Psnc,samp_start);
% datpused1=zeros(Psnc,samp_start);
% datp1=zeros(Psnc,samp_start);
% Istnsnc1=zeros(Psnc,samp_start);
% Isncsnc1=zeros(Psnc,samp_start);
% Istnstn1=zeros(Pstn,samp_start);
% Igpestn1=zeros(Pstn,samp_start);
% Igpegpe1=zeros(Pgpe,samp_start);
% Istngpe1=zeros(Pgpe,samp_start);
% dV_snc1=zeros(Psnc,samp_start);
% dcaer1=zeros(Psnc,samp_start);
% dcamt1=zeros(Psnc,samp_start);
% dcam1=zeros(Psnc,samp_start);
% dcalb1=zeros(Psnc,samp_start);

yeda1T=zeros(1,samp_start);ycda1T=zeros(1,samp_start);
yatpT=zeros(1,samp_start);yrosT=zeros(1,samp_start);
yeda1=zeros(1,samp_start);yatp=zeros(1,samp_start);
yDAT=zeros(1,samp_start);Isp=zeros(1,samp_start);
deda=[];datp=[];dedaT=[];dcdaT=[];datpT=[];drosT=[];
dncS=[];dncT=[];dDAT=[];dIsp=[];
ycaiS=zeros(1,samp_start);ycaiT=zeros(1,samp_start);
ycamtS=zeros(1,samp_start);
dcaiS=[];dcaiT=[];dcamtS=[];
% yrosT=zeros(PsncT,samp_start);
% ycamtS=zeros(Psnc,samp_start);

% yeda1=zeros(1,samp_start);
% dapop1=zeros(Psnc,samp_start);
dnc1=zeros(1,samp_start);
% dDA1=zeros(1,samp_start);
dnc2=zeros(1,samp_start);
% dros_mit1=zeros(Psnc,samp_start);
% dcai1=zeros(Psnc,samp_start);
% datpused1=zeros(Psnc,samp_start);
% datp1=zeros(Psnc,samp_start);
% Istnsnc1=zeros(Psnc,samp_start);
% Isncsnc1=zeros(Psnc,samp_start);
% Istnstn1=zeros(Pstn,samp_start);
% Igpestn1=zeros(1,samp_start);
% Igpegpe1=zeros(1,samp_start);
% Istngpe1=zeros(1,samp_start);
% dV_snc1=zeros(Psnc,samp_start);
% dcaer1=zeros(Psnc,samp_start);
% dcamt1=zeros(Psnc,samp_start);
% dcam1=zeros(Psnc,samp_start);
% dcalb1=zeros(Psnc,samp_start);

NEWstn2snc=1;wDA_snc=1;

%SP Parameters:
tau_f = 200;
tau_r = 10;
% log_e = 0.4342944819;
% e = log_e;
beta = 0.47;
% t = 0;
% lambda_p = 5.5;
Np_t=zeros(Psnc,40/dt);

%%%SELF-KILLING
% idxx = randsample(1:Psnc,Psnc);
% sst=2000/dt;
% ssp=1;
indsapp=[];

%% Energy deficit
enfatra=1.*ones(Msnc,Nsnc);
enfamit=1.*ones(Msnc,Nsnc);
tatp=numel(enfatra);
peratp=round(peren*tatp/100);
idxED = randsample(1:tatp,tatp);
pratp=idxED(1:peratp);
endf=10000/dt;
countt=1;
indsappcamt=[];
indsappcaer=[];

Sig_mts=0*ones(Msnc,Nsnc);
Sig_ers=0*ones(Msnc,Nsnc);

lmin=0.00001;lmax=0.00005;
lam=lmin+rand(Msnc,Nsnc)*(lmax-lmin);

%% Energy deficit terminals
enfatraT=1.*ones(MsncT,NsncT);
enfamitT=1.*ones(MsncT,NsncT);
tatpT=numel(enfatraT);
peratpT=round(perenT.*tatpT./100);
idxEDT = randsample(1:tatpT,tatpT);
pratpT=idxEDT(1:peratpT);
endfT=1000/dt;
counttt=1;
indsappT=[];
lmin=0.00001;lmax=0.00005;
lamT=lmin+rand(MsncT,NsncT)*(lmax-lmin);

% Initiating Therapy after certain percentage cell loss
NcellLoss=round((cl*Psnc)/100); % 50% cell loss
NcellLossT=round((clT*PsncT)/100); % 50% terminal loss

%%%%%% Therrapeutics %%%%%%
% Glutamate Inhibitor
gi=1;
Inc=0;
gi_dose=1;
%         enfatra=0.2*ones(Msnc,Nsnc);%
%         enfamit
% clit=cell(1,Ttime);
limsg=0.1;
limsm=0.15;
limsgT=0.1;
limsmT=0.15;
k2=1;
% Imsn1_snc=0;Imsn2_snc=0;
con=1;
% Isp=[];
bre=0;
ldopaS=0;
ldopaT=0;
SPT=1;
% sLD=0;
% sLDT=0;
%%
forstart1=forstart;
for k = 1:Ttime
    
    %     if k>k2
    %         if Inc>=NcellLoss
    %             gi=gi_dose;
    %         end
    %     end
    %     %LDOPA medication
    if LDOPAon==1
        if k>k2
            if Inc>=NcellLoss && Son==1
                sLD=ldopaS_dose;
                sLDT=ldopaT_dose;
            end
            
            if IncT>=NcellLossT && Ton==1
                sLD=ldopaS_dose;
                sLDT=ldopaT_dose;
            end
        end
    end
    
    if SPon==1
        if k>k2
            if Inc>=NcellLoss
                SPT=spt_dose;
            end
            
            if IncT>=NcellLossT
                SPT=spt_dose;
            end
        end
    end
    
    if GSon==1
        if k>k2
            if Inc>=NcellLoss
                GSH_totT=gst_dose;
                GSHT=gst_dose;
            end
            
            if IncT>=NcellLossT
                GSH_totT=gst_dose;
                GSHT=gst_dose;
            end
        end
    end
    
    %     ATP=2.4*ones(Msnc,Nsnc);
    
    %     if k==sst
    %         indsapp=idxx(1:ssp);
    %         ssp=ssp+1;
    %         sst=sst+(200/dt);
    %     end
    
    V_snc(indsapp) = -80.*ones(size(indsapp));
    apop(indsapp)=zeros(size(indsapp));
    Ca_i(indsapp)=zeros(size(indsapp)); %mM
    eda(indsapp) = 26e-6.*zeros(size(indsapp)); %mM
    ATPused(indsapp)=0.*zeros(size(indsapp));
    ATP(indsapp)=0.*zeros(size(indsapp));
    
    if k>endf
        Renfatra=1.*exp(-countt.*lam);
        Renfamit=1.*exp(-countt.*lam);
        %         edda=1e-5.*exp(-counttt.*lam);
        countt=countt+1;
        
        enfatra(pratp)=Renfatra(pratp);
        enfamit(pratp)=Renfamit(pratp);
        %         eda(pratp)=edda(pratp);
        
        enfatra(enfatra<limsg)=limsg;
        enfamit(enfamit<limsm)=limsm;
        %        nATP(nATP<rATP1)=rATP1;
        
    end
    
    %     if(k > 5000/dt)
    %         enfatra=0.2*ones(Msnc,Nsnc);%
    %         enfamit=0.1*ones(Msnc,Nsnc);;
    %     else
    %         enfatra=1*ones(Msnc,Nsnc);;
    %         enfamit=1*ones(Msnc,Nsnc);;
    %     end
    %
    
    %     phier=Ca_i-Ca_er;
    
    %     inds1=Ca_mt(pratp)>camtthr;
    %     inds2=pratp(inds1);
    
    inds2=find(Ca_mt>camtthr);
    indsappcamt=[indsappcamt inds2];
    if isempty(inds2)==0
        indsappcamt=unique(indsappcamt);
        inds2=[];
    end
    Ca_mt(indsappcamt)=zeros(size(indsappcamt)); %mM
    Ca_mt(isnan(Ca_mt)| isinf(Ca_mt)| Ca_mt<0)=0;
    %
    %     indsapcaer=find(phier>0.0);
    %     indsappcaer=[indsappcaer indsapcaer'];
    %     indsappcaer=unique(indsappcaer);
    %
    Sig_mts(indsappcamt)=0.01.*ones(size(indsappcamt));
    %     Sig_ers(indsappcaer)=0.01.*ones(size(indsappcaer));
    
    DA = RescaleRange(eda1,1e-5,2.1e-5,0,1);
    
    DA(isnan(DA)| isinf(DA)| DA<0)=0;
    
    if DA<0
        DA=0;
    end
    if DA>1
        DA=1;
    end
    
    
    wda_gpe=1;
    %     ssmax = RescaleRange(DA,0,1,40,0.1);
    wsg=((1-CD2*da))*wsg1;
    wgs=((1-CD2*DA))*wgs1;
    wlatstn = weightcal_stn(DA*pd);
    V_snc(isnan(V_snc)| isinf(V_snc))=0;
    Ssnc(isnan(Ssnc)| isinf(Ssnc))=0;
    
    %----------------------------------------SNc-----------------------------------------%
    % Lateral SNc-SNc connections
    wlatsnc = weightcal_snc(DA,rss,nlat);
    Hsnc=1./(1+exp(-(V_snc-20+57)./2));
    Ssncnxt=Ssnc+((2.*Hsnc.*(1-Ssnc))-Ssnc.*0.08).*dt;
    Ssncfin=conv2(Ssncnxt,wlatsnc,'same');
    Igabasnc=1.*0.01.*(V_snc+63.45).*Ssncfin;
    Igabasnc(isnan(Igabasnc)| isinf(Igabasnc))=0;
    
    % STN-SNc connections
    if (NEWstn2snc==1)
        start=1;stop=nproj;totspk=zeros(1,Psnc);
        for ll=1:Psnc
            %             stn_snccurr1=zeros(proj,proj);
            %             stn_snccurr1=reshape(stncurr_spk(idx(start:stop)),proj,proj);
            totspk(ll)=sum(stncurr_spk(idx(start:stop)));
            start=start+nproj;stop=stop+nproj;
        end
    end
    stn_snccurr=zeros(8,8);
    stn_snccurr=reshape(totspk,8,8);
    
    h_nmdasnc = (1-lam_nmda).* h_nmdasnc + lam_nmda.*stn_snccurr; %psp nmda snc
    h_ampasnc = (1-lam_ampa).* h_ampasnc + lam_ampa.*stn_snccurr;
    
    tmp_nmda_snc = h_nmdasnc.*(Enmda - V_snc);
    tmp_ampa_snc = h_ampasnc.*(Eampa - V_snc);
    
    
    I_nmda_snc = wDA_snc.*gi.*wstnsnc_matrix.*tmp_nmda_snc;
    I_ampa_snc = wDA_snc.*gi.*wstnsnc_matrix.*tmp_ampa_snc;
    
    B = 1./(1 + (mg0./3.57).*exp(-0.062.*V_snc));
    
    Istnsnc=1*(B.*I_nmda_snc + I_ampa_snc);
    Istnsnc(isnan(Istnsnc)| isinf(Istnsnc))=0;
    
    
    %MSN1 to SNc------%
    
    spkmsn1=sum(msn1curr_spk(idx_MSN_SNc(:,:)),2);
    spkmsn11=reshape(spkmsn1,Msnc,Nsnc);
    
    hlat_msn1_snc = (1-lam_gaba).* hlat_msn1_snc + lam_gaba.*spkmsn11;
    Imsn1_snc = w_msn1_snc.*(hlat_msn1_snc.*(V_snc - Egaba));
    %     I_msn1_snc=[I_msn1_snc Imsn1_snc];
    Imsn1_snc(isnan(Imsn1_snc)| isinf(Imsn1_snc))=0;
    
    %MSN2 to SNc------%
    
    spkmsn2=sum(msn2curr_spk(idx_MSN_SNc(:,:)),2);
    spkmsn22=reshape(spkmsn2,Msnc,Nsnc);
    
    hlat_msn2_snc = (1-lam_gaba).* hlat_msn2_snc + lam_gaba.*spkmsn22;
    Imsn2_snc = w_msn2_snc.*(hlat_msn2_snc.*(V_snc - Egaba));
    %     I_msn2_snc=[I_msn2_snc Imsn2_snc];
    Imsn2_snc(isnan(Imsn2_snc)| isinf(Imsn2_snc))=0;
    
    % SP-neuropeptide
    %MSN2------%
    var1 = -spkmsn2./tau_f;
    var2 = -spkmsn2./tau_r;
    var3 = exp(var1);
    var4 = exp(var2);
    Ap_t = (var3 - var4);
    
    var_5 = ((-Ap_t./5.5));
    Np_t0 = beta.*(1 - (exp(var_5)).^(2.5));
    
    if k<=40/dt
        Np_t(:,k)=Np_t0;
        Isp_snc =Istnsnc;
    else
        Np_t(:,1)=[];
        Np_t(:,k-con)=Np_t0;
        Nsp=reshape(Np_t(:,1),Msnc,Nsnc);
        DASP=(1-DA_SP.*mean2(phi1));
        DASP(DASP<0)=0.00001;
        Isp_snc = Istnsnc.*DASP.*(1+wsp.*SPT.*Nsp);
        con=con+1;
    end
    Isp_snc(isnan(Isp_snc)| isinf(Isp_snc))=0;
    
    
    
    %     Isp=[Isp sum(sum(Isp_snc))];
    
    %     Iz_ms = 0.00003.*(1 + 1.*(Np_t.*(tspan(k)-40)));
    %     if k>40/dt
    %         Np_t(:,1)=[];
    %         Nsp=reshape(Np_t(:,1),Msnc,Nsnc);
    %         Isp_msn = Istnsnc.*500.*Nsp;
    % %         Isp_snc = Istnsnc.*(1 + (1-0.7.*noSC).*(1-DA_SP.*phi1).*wsp.*10000.*Np_t(k-40/dt));
    %     else
    %         Isp_snc =Istnsnc;
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Stimulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Ibg=0;
    if(k < del + dur && k > del)
        Iapp = Ibg +Istim;
    else
        Iapp = Ibg;
    end
    Iext=Iapp;
    
    %Soma
    
    %ATP-dependent equations
    k_1pc = 1.00000./(1.00000+0.100000./ATP);
    k_1nk = 0.370000./(1.00000+0.0940000./ATP);
    
    %Soma
    %Membrane potential
    VD = V_snc./V_tau;
    
    % HCN current
    kf_free = 0.00600000./(1.00000+exp((V_snc+87.7000)./6.45000));
    kf_bnd = 0.0268000./(1.00000+exp((V_snc+94.2000)./13.3000));
    kf_hcn = kf_free.*P_c+kf_bnd.*(1.00000-P_c);
    kr_free = 0.0800000./(1.00000+exp(-(V_snc+51.7000)./7.00000));
    kr_bnd = 0.0800000./(1.00000+exp(-(V_snc+35.5000)./7.00000));
    kr_hcn = kr_free.*P_o+kr_bnd.*(1.00000-P_o);
    
    % Calcium binding proteins
    CaCalb = Calbtot-Calb;
    J_calb = kcal_1.*Calb.*Ca_i-kcal_2.*CaCalb;
    
    CaCam = Camtot-Cam;
    kcam_cb = 12000.0.*(power(Ca_i, 2.00000));
    kcam_nb = 3.70000e+06.*(power(Ca_i, 2.00000));
    alpha_cam = kcam_cb.*kcam_nb.*(1.00000./(kcam_cb+kcam_nd)+1.00000./(kcam_cd+kcam_nd));
    beta_cam = kcam_cd.*kcam_nd.*(1.00000./(kcam_cb+kcam_nd)+1.00000./(kcam_cd+kcam_nd));
    J_cam = alpha_cam.*Cam-beta_cam.*CaCam;
    
    
    K_pci = (173.600./(1.00000+CaCam./5.00000e-05)+6.40000).*1.00000e-05;
    P_E1Spc = 1.00000./(1.00000+K_pci./Ca_i);
    P_E1pc = 1.00000-P_E1Spc;
    alpha_pc = k_1pc.*P_E1Spc+k_3pc.*P_E1pc;
    
    V_Ca = 0.500000.*log(Ca_o./Ca_i);
    h_cal = 0.000450000./(0.000450000+Ca_i);
    I_CaL = (g_cal.*m_cal.*h_cal.*(power(Ca_i.*Ca_o, 1.0./2)).*sinh(VD-V_Ca))./(sinh(VD)./VD);
    K_pmca = k_pmca.*((10.5600.*CaCam)./(CaCam+5.00000e-05)+1.20000);
    I_pmca = K_pmca.*(k_1pc.*P_E1Spc.*y_pc-k_2pc.*P_E2Spc.*(1.00000-y_pc)).*1.00000;
    Dr = (1.00000+0.00100000.*((power(Na_i, 3.00000)).*Ca_o+(power(Na_o, 3.00000)).*Ca_i)).*(1.00000+Ca_i./0.00690000);
    I_xm = (k_xm.*((power(Na_i, 3.00000)).*Ca_o.*exp(dell.*VD)-(power(Na_o, 3.00000)).*Ca_i.*exp((dell-1.00000).*VD)))./Dr;
    J_ca = (-1.00000./(2.00000.*F.*vol_cyt)).*((I_CaL+2.00000.*I_pmca)-2.00000.*I_xm);
    
    V_Na = log(Na_o./Na_i);
    O_na = (power(m_na, 3.00000)).*h_na;
    I_Na = (g_na.*O_na.*(power(Na_i.*Na_o, 1.0./2)).*sinh(0.500000.*(VD-V_Na)))./(sinh(0.500000.*VD)./(0.500000.*VD));
    I_Nalk = (g_nalk.*(power(Na_i.*Na_o, 1.0./2)).*sinh(0.500000.*(VD-V_Na)))./(sinh(0.500000.*VD)./(0.500000.*VD));
    I_NaHCN = (g_nahcn.*O_hcn.*(power(Na_i.*Na_o, 1.0./2)).*sinh(0.500000.*(VD-V_Na)))./(sinh(0.500000.*VD)./(0.500000.*VD));
    P_E1Snk = 1.00000./(1.00000+(K_nknai./Na_i).*(1.00000+K_i./K_nkki));
    Na_eff = Na_o.*exp(-0.820000.*VD);
    P_E2Snk = 1.00000./(1.00000+(K_nknao./Na_eff).*(1.00000+K_o./K_nkko));
    I_nk = k_nk.*(k_1nk.*P_E1Snk.*y_nk-k_2nk.*P_E2Snk.*(1.00000-y_nk)).*1.00000;
    J_Na = (-1.00000./(F.*vol_cyt)).*(3.00000.*I_nk+3.00000.*I_xm+I_Na+I_Nalk+I_NaHCN);
    
    P_E1Dnk = 1.00000./(1.00000+(K_nkki./K_i).*(1.00000+Na_i./K_nknai));
    alpha_nk = k_1nk.*P_E1Snk+k_3nk.*P_E1Dnk;
    P_E2Dnk = 1.00000./(1.00000+(K_nkko./K_o).*(1.00000+Na_eff./K_nknao));
    beta_nk = k_2nk.*P_E2Snk+k_4nk.*P_E2Dnk;
    
    V_K = log(K_o./K_i);
    O_sk = (power(Ca_i, 4.20000))./(power(0.000350000, 4.20000)+power(Ca_i, 4.20000));
    I_Ksk = (g_ksk.*O_sk.*(power(K_i.*K_o, 1.0./2)).*sinh(0.500000.*(VD-V_K)))./(sinh(0.500000.*VD)./(0.500000.*VD));
    O_kdr = power(m_kdr, 3.00000);
    I_Kdr = g_kdr.*O_kdr.*(V_snc-V_K.*V_tau);
    O_kir = 1.00000./(1.00000+exp((V_snc+85.0000)./12.1000));
    I_Kir = g_kir.*O_kir.*(V_snc-V_K.*V_tau);
    I_K = I_Ksk+I_Kdr+I_Kir;
    J_K = (-1.00000./(F.*vol_cyt)).*(I_K-2.00000.*I_nk);
    
    % ER
    J_pump = k_pump.*Ca_i.*ATP; %J_pump
    J_ch = k_ch.*((power(Ca_i, 2.00000))./(power(K1, 2.00000)+power(Ca_i, 2.00000))).*(Ca_er-Ca_i); %J_ch
    J_leak = k_leak.*(Ca_er-Ca_i); %J_leak
    
    % Mito
    J_out = (k_out.*((power(Ca_i, 2.00000))./(power(K3, 2.00000)+power(Ca_i, 2.00000)))+k_m).*Ca_mt; % J_out
    J_in = k_in.*((power(Ca_i, 8.00000))./(power(K2, 8.00000)+power(Ca_i, 8.00000)));%.*ATP; % J_in
    
    % Calcium dynamics
    J_Ca = J_ca-(J_calb+4.00000.*J_cam)-J_pump+J_ch+J_leak-J_in+J_out;
    
    adca=0;
    %Terminal
    
    % ATP-dependent DA packing
    %     ada=RescaleRange(ATP,0.2,2.3,0.001,1);
    ada=0.001.*(exp(3.*ATP));
    
    % ATP-dependent vescile recycling
    %     nRRP=5;
    eta_nrrp=eta_nrrp_max-beta_nrrp_asyn.*(1./(1+(power((Kasyn./ASYNA),4))));
    eta_nrrp(isnan(eta_nrrp)| isinf(eta_nrrp)| eta_nrrp<0)=0;
    nRRP=1.*eta_nrrp.*(exp(1.*ATP));
    
    Ca_i(isnan(Ca_i)| isinf(Ca_i)| Ca_i<0)=0;
    
    Vsynt = Vsynt_max./(((Ksynt./(adca+Ca_i)).^4)+1);
    jsynt = (Vsynt./(1+((Ktyr./TYR).*(1+(cda./Kicda)+(eda./Kieda)))));
    
    jvmat = ada.*Vcda_max .* MM_kin(cda,Kcda,1);%.*ATP;
    
    jida = kmao .* cda;
    
    %     nRRP = (40./((1+exp(-(vda-vdao)./vdas)).*(1+exp((eda-dara)./dars))));
    prob = 0.14 .* MM_kin((adca+Ca_i),krel,4);
    jrel = psi .* nRRP .* prob;
    
    jdat = Veda_max .* MM_kin(eda,Keda,1);
    
    jeda = kcomt .* eda;
    
    jldopa = Vaadc_max .* MM_kin(LDOPA,Kaadc,1);
    
    % Energy metabolism
    % Energy consumed by active pumps
    % V_pumps=0;
    V_pumps1 = 1.*(1.00000./(F.*vol_cyt)).*(I_nk+I_pmca);
    V_pumps2 = 1.*(jvmat);
    V_pumps3 = 100.*jrel;
    
    v_stim=0;
    % v_stim1=(1.00000/(F*vol_cyt))*(I_nk+I_pmca);
    v_stim1=0.0.*(I_nk+I_pmca);
    J_er=(beta_er./rho_er).*(J_pump);
    
    %     V_id(k)=V_pumps1;
    %     V_dp(k)=V_pumps2;
    %     V_er(k)=J_er;
    
    V_pumps=V_pumps1+V_pumps2+J_er+V_pumps3;
    
    uADP = power(Q_adk, 2.00000)+4.00000.*Q_adk.*(ANP./ATP-1.00000);
    ADP = (ATP./2.00000).*(-Q_adk+power(uADP, 1.0./2));
    Cr = PCr_tot-PCr;
    V_ck = 0.*(kf_ck.*PCr.*ADP-kr_ck.*Cr.*ATP);
    
    ATP_inh = power((1.00000+nATP.*(ATP./KI_ATP))./(1.00000+ATP./KI_ATP), 4.00000);
    ATP_inh(isnan(ATP_inh)| isinf(ATP_inh)| ATP_inh<0)=0;
    V_pk = enfatra.*Vmax_pk.*(GAP./(GAP+Km_GAP_pk)).*(ADP./(ADP+Km_ADP_pk)).*ATP_inh;
    pa=1;
    V_op = enfamit.*Vmax_op.*((pa.*PYR)./((pa.*PYR)+Km_PYR_op)).*(ADP./(ADP+Km_ADP_op)).*(1.00000./(1.00000+0.100000.*(ATP./ADP)));
    
    AMP = ANP-(ATP+ADP);
    AMP_act = power((1.00000+AMP./Ka_AMP)./(1.00000+nAMP.*(AMP./Ka_AMP)), 4.00000);
    V_pfk = Vmax_pfk.*(F6P./(F6P+Km_F6P_pfk)).*(ATP./(ATP+Km_ATP_pfk)).*(F26P./(F26P+Km_F26P_pfk)).*ATP_inh.*AMP_act;
    
    eta_ldh=1-beta_ldh_ros.*((ROS.^4)./((ROS.^4)+(Kldh_ros.^4)));
    eta_ldh(isnan(eta_ldh)| isinf(eta_ldh)| eta_ldh<0)=0;
    V_ldh = 1.*eta_ldh.*(kf_ldh.*PYR-kr_ldh.*LAC);
    V_lac = Vlac_0.*(1.00000+v_stim1.*K_lac)-K_lac_eff.*LAC;
    
    V_hk = Vmax_hk.*(ATP./(ATP+Km_ATP_hk)).*(power(1.00000+power(F6P./KI_F6P, 4.00000), -1.00000)).*GLCe;
    AMP_pfk2 = (power(AMP./Kamp_pfk2, nh_amp))./(1.00000+power(AMP./Kamp_pfk2, nh_amp));
    V_pfk2 = Vmaxf_pfk2.*(ATP./(ATP+Km_ATP_pfk2)).*(F6P./(F6P+Km_F6P_pfk2)).*AMP_pfk2-Vmaxr_pfk2.*(F26P./(F26P+Km_F26P_pfk2));
    
    V_ATPase = Vmax_ATPase.*(ATP./(ATP+Km_ATP)).*(1.00000+v_stim);
    dAMP_dATP = -1.00000+Q_adk./2.00000+-(0.500000.*(power(uADP, 1.0./2)))+Q_adk.*(ANP./(ATP.*(power(uADP, 1.0./2))));
    dAMP_dATP(isnan(dAMP_dATP)| isinf(dAMP_dATP)| dAMP_dATP<0)=0;
    
    eta_op=eta_op_max-beta_op_asyn.*(((ASYNA.^4)./((ASYNA.^4)+(Kasyn.^4))));
    eta_op(isnan(eta_op)| isinf(eta_op)| eta_op<0)=0;
    
    GSSG=(GSH_tot-GSH)./2;
    NADP=NADPH_tot-NADPH;
    Vppp = Vmax_ppp.*(F6P./(F6P+Km_F6P_pfk))./(1+((NADPH./NADP)./Ki_nadph));
    Vgr = kf_gr.*GSSG.*NADPH-kr_gr.*GSH.*NADP;
    
    Vros_leak=(0.5282./ATP).*(1-eta_op).*V_op;%0.221 or (0.5282./ATP)
    %     Vros_leak=0.221.*(1-eta_op).*V_op;%0.221 or (0.5282./ATP)
    
    Vros_cat=Kros_cat.*ROS;
    
    %     Vros_dopa=0.*Kros_dopa.*(1./(1+power(Kasyn./ASYNA,4)));%(((ASYNA.^4)./((ASYNA.^4)+(Kasyn.^4))));
    Vros_dopa=0.001.*Kros_dopa.*(((cda)./((cda)+(Kasyn.*1000))));
    
    Vros_dox=Kros_dox.*GSH.*ROS;
    
    Vasyn_syn=Kasyn_syn;
    
    Vasyn_ox=Kasyn_ox.*ROS.*ASYN;
    
    Vasyn_to=Kasyn_to.*ASYN;
    
    Vasyn_agg=Krasyn_agg.*ASYNA.*(1./(1+power(Kasyn_agg./ASYNA,6)));%(((ASYNA.^6)./((ASYNA.^6)+(Kasyn_agg.^6))));
    
    Ub=Ub_tot-ASYNT;
    Vasyn_tag=Kasyn_tag.*ASYNA.*Ub.*ATP;
    
    Vasyn_prt=Krasyn_prt.*ASYNT.*ATP.*(1-beta_asyn_prt.*(1./(1+power(Kasyn_prt./ASYNG,4))));%(((ASYNG.^4)./((ASYNG.^4)+(Kasyn_prt.^4)))));
    
    Vasyn_lyso=Kasyn_lyso.*ASYNG.*ATP;
    
    Vasyn_lb=Krasyn_lb.*ASYNG.*(1./(1+power(Kasyn_lb./ASYNG,6)));%(((ASYNG.^6)./((ASYNG.^6)+(Kasyn_lb.^6))));
    
    vAADCext = 0.24.*(vAADC_Vmax)./(1+(power((vAADC_Km./sLD),0.5)));
    
    ext=(power(1-dAMP_dATP, -1.00000));
    ext((isnan(ext)| isinf(ext)| ext<0))=0;
    
    %%%%%%%%%%%%%%%%%%%%%% Differential equations %%%%%%%%%%%%%%%%%%%%%%%%%
    
    V_sncnxt = V_snc + (((F.*vol_cyt)./(C_sp.*A_pmu)).*(J_Na+J_K+2.00000.*J_Ca+Iext-Igabasnc+(scfa.*Isp_snc)-(scfamsn.*Imsn1_snc)-(scfamsn.*Imsn2_snc))).*dt;
    Ca_inxt = Ca_i + (J_Ca).*dt;
    Na_inxt = Na_i + (J_Na).*dt;
    K_inxt = K_i + (J_K).*dt;
    Calbnxt = Calb + (-J_calb).*dt;
    Camnxt = Cam + (-J_cam).*dt;
    m_calnxt = m_cal + ((1.00000./(1.00000+exp(-(V_snc+15.0000)./7.00000))-m_cal)./(7.68000.*exp(-(power((V_snc+65.0000)./17.3300, 2.00000)))+0.723100)).*dt;
    m_nanxt = m_na + (A_mna.*exp(za_mna.*VD).*(1.00000-m_na)-B_mna.*exp(-zb_mna.*VD).*m_na).*dt;
    h_nanxt = h_na + (A_hna.*exp(za_hna.*VD).*(1.00000-h_na)-B_hna.*exp(-zb_hna.*VD).*h_na).*dt;
    O_hcnnxt = O_hcn + (kf_hcn.*(1.00000-O_hcn)-kr_hcn.*O_hcn).*dt;
    m_kdrnxt = m_kdr + ((1.00000./(1.00000+exp(-(V_snc+25.0000)./12.0000))-m_kdr)./(18.0000./(1.00000+exp((V_snc+39.0000)./8.00000))+1.00000)).*dt;
    y_pcnxt = y_pc + (beta_pc.*(1.00000-y_pc)-alpha_pc.*y_pc).*dt;
    y_nknxt = y_nk + (beta_nk.*(1.00000-y_nk)-alpha_nk.*y_nk).*dt;
    ATPusednxt = ATPused+(-ATPused+(1.00000./(F.*vol_cyt)).*(I_nk+I_pmca)).*dt;%ATPused
    Ca_ernxt = Ca_er + ((beta_er./rho_er).*(J_pump-(J_ch+J_leak))).*dt;
    Ca_mtnxt = Ca_mt + ((beta_mt./rho_mt).*(J_in-J_out)).*dt;
    cdanxt = cda+(0*jsynt + jdat - jvmat - jida + jldopa).*dt;%cda
    vdanxt = vda+(jvmat - jrel).*dt;%vda
    edanxt = eda+(jrel - jdat - jeda + vAADCext).*dt;%eda
    calnxt = cal+(-k3f.*(Sig_ers.*cal)+k3b.*(cai_cal)).*dt;%cal
    cai_calnxt = cai_cal+(k3f.*(Sig_ers.*cal)-k3b.*(cai_cal)-k4f.*(cai_cal)).*dt;%cai_cal
    cal_actnxt = cal_act+(k4f.*(cai_cal)-k5f.*(cal_act.*casp12)+k5b.*(cal_act_casp12)).*dt;%cal_act
    casp12nxt = casp12+(-k5f.*(cal_act.*casp12)+k5b.*(cal_act_casp12)).*dt;%casp12
    cal_act_casp12nxt = cal_act_casp12+(k5f.*(cal_act.*casp12)-k5b.*(cal_act_casp12)-k6f.*(cal_act_casp12)).*dt;%cal_act_casp12
    casp12_actnxt = casp12_act+(k6f.*(cal_act_casp12)-k7f.*(casp12_act.*casp9)+k7b.*(casp12_act_casp9)).*dt;%casp12_act
    casp9nxt = casp9+(-k7f.*(casp12_act.*casp9)+k7b.*(casp12_act_casp9)).*dt;%casp9
    casp12_act_casp9nxt = casp12_act_casp9+(k7f.*(casp12_act.*casp9)-k7b.*(casp12_act_casp9)-k8f.*(casp12_act_casp9)).*dt;%casp12_act_casp9
    casp9_actnxt = casp9_act+(k8f.*(1.*casp12_act_casp9)+k9b.*(casp9_act_casp3)-k9f.*(casp9_act.*casp3)+1.*k28f.*Cytc_casp9-k12f.*casp9_act.*IAP+k12b.*casp9_act_IAP).*dt;%casp9_act
    casp3nxt = casp3+(-k9f.*(casp9_act.*casp3)+k9b.*(casp9_act_casp3)).*dt;%casp3
    casp9_act_casp3nxt = casp9_act_casp3+(-k10f.*(casp9_act_casp3)-k9b.*(casp9_act_casp3)+k9f.*(casp9_act.*casp3)).*dt;%casp9_act_casp3
    casp3_actnxt = casp3_act+(k10f.*(casp9_act_casp3)-k11f.*(casp9_act.*casp3_act)-k13f.*casp3_act.*IAP+k13b.*casp3_act_IAP).*dt;%casp3_act
    apopnxt =  apop+(k11f.*(casp9_act.*casp3_act)).*dt;%apop
    
    F6Pnxt = F6P+(V_hk-(V_pfk-V_pfk2)).*dt;%F6P
    F26Pnxt = F26P+(V_pfk2).*dt;%F26P
    GAPnxt = GAP+(2.00000.*V_pfk-V_pk).*dt;%GAP
    PYRnxt = PYR+(V_pk-(V_op+V_ldh)).*dt;%PYR
    LACnxt = LAC+(2.25000.*V_ldh+V_lac).*dt;%LAC
    ATPnxt = ATP+(((1.*(1.*2.00000.*V_pk+15.0000.*eta_op.*V_op+V_ck))-(V_hk+V_pfk+V_pfk2+V_ATPase+V_pumps)).*ext).*dt;%ATP
    PCrnxt = PCr+(-V_ck).*dt;%PCr
    
    ROS_mitnxt = ROS_mit+(k29f.*Sig_mts.*Mit).*dt;%ROS_mit
    PTP_mit_actnxt = PTP_mit_act+(k30f.*ROS_mit.*PTP_mit).*dt;%PTP_mit_act
    Cytc_mitnxt = Cytc_mit+(-k31f.*PTP_mit_act.*Cytc_mit).*dt;%Cytc_mit
    Cytcnxt = Cytc+(-k27f.*Cytc.*casp9+k27b.*Cytc_casp9+k31f.*PTP_mit_act.*Cytc_mit).*dt;%Cytc
    Cytc_casp9nxt = Cytc_casp9+(k27f.*Cytc.*casp9-k27b.*Cytc_casp9-k28f.*Cytc_casp9).*dt;%Cytc_casp9
    IAPnxt = IAP+(-k12f.*casp9_act.*IAP+k12b.*casp9_act_IAP-k13f.*casp3_act.*IAP+k13b.*casp3_act_IAP).*dt;%IAP
    casp9_act_IAPnxt = casp9_act_IAP+(k12f.*casp9_act.*IAP-k12b.*casp9_act_IAP).*dt;%casp9_act_IAP
    casp3_act_IAPnxt = casp3_act_IAP+(k13f.*casp3_act.*IAP-k13b.*casp3_act_IAP).*dt;%casp3_act_IAP
    
    ROSnxt = ROS+(Vros_leak + Vros_ex - Vros_cat + Vros_dopa - Vros_dox).*dt;%ROS
    ASYNnxt = ASYN+(Vasyn_syn - Vasyn_ox - Vasyn_to).*dt;%ASYN
    ASYNAnxt = ASYNA+(Vasyn_ox - Vasyn_agg - Vasyn_tag).*dt;%ASYNA
    ASYNTnxt = ASYNT+(Vasyn_tag - Vasyn_prt).*dt;%ASYNT
    ASYNGnxt = ASYNG+(Vasyn_agg - Vasyn_lyso - Vasyn_lb).*dt;%ASYNG
    LBnxt = LB+(Vasyn_lb).*dt;%LB
    
    NADPHnxt = NADPH+(2.*Vppp - Vgr).*dt;%NADPH
    GSHnxt = GSH+(2.*Vgr - 2.*Vros_dox).*dt;%
    
    LDOPAnxt = LDOPA+(((Vtran_max.*sLD)./(Ksld.*(1+(sTYR./Kstyr)+(sTRP./Kstrp))+sLD)) - jldopa + jsynt).*dt;%ldopa
    
    
    %     inds=find(V_snc <=80 & V_snc >=-20);
    %     snc_firings=[snc_firings; k+0*inds,inds+0*inds];
    
    %     if k>sncstart/dt
    %         [snc_firings]=ConvertAPtoST(snc_firings,Psnc);
    %         sncstart=sncstart+500;
    %     end
    
    V_snc=V_sncnxt;m_cal=m_calnxt;m_kdr=m_kdrnxt;m_na=m_nanxt;
    h_na=h_nanxt;O_hcn=O_hcnnxt;Calb=Calbnxt;
    Cam=Camnxt;y_nk=y_nknxt;y_pc=y_pcnxt;
    K_i=K_inxt;Na_i=Na_inxt;Ca_i=Ca_inxt;
    Ca_er=Ca_ernxt;Ca_mt=Ca_mtnxt;
    cda=cdanxt;vda=vdanxt;eda=edanxt;ATPused=ATPusednxt;
    cal=calnxt;cai_cal=cai_calnxt;cal_act=cal_actnxt;
    casp12=casp12nxt;cal_act_casp12=cal_act_casp12nxt;casp12_act=casp12_actnxt;
    casp9=casp9nxt;casp12_act_casp9=casp12_act_casp9nxt;casp9_act=casp9_actnxt;
    casp3=casp3nxt;casp9_act_casp3=casp9_act_casp3nxt;casp3_act=casp3_actnxt;
    apop=apopnxt;
    ROS_mit=ROS_mitnxt;
    PTP_mit_act=PTP_mit_actnxt;
    Cytc_mit=Cytc_mitnxt;
    Cytc=Cytcnxt;
    Cytc_casp9=Cytc_casp9nxt;
    IAP=IAPnxt;
    casp9_act_IAP=casp9_act_IAPnxt;
    casp3_act_IAP=casp3_act_IAPnxt;
    F6P=F6Pnxt;
    F26P=F26Pnxt;
    GAP=GAPnxt;
    PYR=PYRnxt;
    LAC=LACnxt;
    ATP=ATPnxt;
    PCr=PCrnxt;
    Ssnc=Ssncnxt;
    ROS=ROSnxt;
    ASYN=ASYNnxt;
    ASYNA=ASYNAnxt;
    ASYNT=ASYNTnxt;
    ASYNG=ASYNGnxt;
    LB=LBnxt;
    LDOPA=LDOPAnxt;
    NADPH=NADPHnxt;
    GSH=GSHnxt;
    
    Ca_i(isnan(Ca_i)| isinf(Ca_i)| Ca_i<0)=0;
    K_i(isnan(K_i)| isinf(K_i)| K_i<0)=0;
    Na_i(isnan(Na_i)| isinf(Na_i)| Na_i<0)=0;
    
    
    %% SNc soma to terminal
    %----------------------------Soma to Terminal-------------------------%
    start=1;stop=nproj;totspk=zeros(1,PsncT);
    for i=1:numel(Ca_i)
        totspk(idx_terminal(start:stop))=Ca_i(i);
        start=start+nproj;stop=stop+nproj;
    end
    %     Ca_iT=zeros(MsncT,NsncT);
    Ca_iT=reshape(totspk,MsncT,NsncT);
    Ca_iT(isnan(Ca_iT)| isinf(Ca_iT)| Ca_iT<0)=0;
    
    Ca_iT(indsappT)=1.*zeros(size(indsappT));
    cdaT(indsappT) = 1e-4.*zeros(size(indsappT)); %mM
    vdaT(indsappT) = 500.*zeros(size(indsappT)); %mM
    edaT(indsappT) = 26e-6.*zeros(size(indsappT)); %mM
    ATPT(indsappT)=0.*zeros(size(indsappT));
    
    %% SNc terminal
    %-----------------------------SNc Terminal----------------------------%
    
    if k>endfT
        Renfatra=1.*exp(-counttt.*lamT);
        Renfamit=1.*exp(-counttt.*lamT);
        counttt=counttt+1;
        
        enfatraT(pratpT)=Renfatra(pratpT);
        enfamitT(pratpT)=Renfamit(pratpT);
        
        enfatraT(enfatraT<limsgT)=limsgT;
        enfamitT(enfamitT<limsmT)=limsmT;
    end
    %     enfatraT=gl;
    %     enfamitT=mt;
    
    adca=0;
    %Terminal
    VsyntT = Vsynt_max./(power((Ksynt./(adca+Ca_iT)),4)+1);
    jsyntT = (VsyntT./(1+((Ktyr./TYR).*(1+(cdaT./Kicda)+(edaT./Kieda)))));
    
    jvmatT = adaT.*Vcda_max .* MM_kin(cdaT,Kcda,1);%.*(ATP);
    
    jidaT = kmao .* cdaT;
    
    %     nRRP = (40./((1+exp(-(vda-vdao)./vdas)).*(1+exp((eda-dara)./dars))));
    %     nRRP = 1;
    
    % ATP-dependent DA packing
    %     ada=RescaleRange(ATP,0.2,2.3,0.001,1);
    adaT=0.001.*(exp(3.*ATPT));
    
    % ATP-dependent vescile recycling
    %     nRRP=10;
    eta_nrrpT=eta_nrrp_max-beta_nrrp_asyn.*(1./(1+(power((Kasyn./ASYNAT),4))));
    eta_nrrpT(isnan(eta_nrrpT)| isinf(eta_nrrpT)| eta_nrrpT<0)=0;
    nRRPT=1.*eta_nrrpT.*(exp(1.*ATPT));
    
    probT = 0.14 .* MM_kin(Ca_iT,krel,4);
    jrelT = psi .* nRRPT .* probT;
    
    jdatT = Veda_max .* MM_kin(edaT,Keda,1);
    
    jedaT = kcomt .* edaT;
    
    jldopaT = Vaadc_max .* MM_kin(LDOPAT,Kaadc,1);
    
    % Energy metabolism
    % Energy consumed by active pumps
    % V_pumps=0;
    V_pumps2T = 1.*(jvmatT);
    V_pumps3T = 100.*jrelT;
    v_stim=0;
    % v_stim1=(1.00000./(F.*vol_cyt)).*(I_nk+I_pmca);
    v_stim1=0;
    
    
    V_pumpsT=V_pumps2T+V_pumps3T;
    
    uADPT = power(Q_adk, 2.00000)+4.00000.*Q_adk.*(ANP./ATPT-1.00000);
    ADPT = (ATPT./2.00000).*(-Q_adk+power(uADPT, 1/2));
    CrT = PCr_tot-PCrT;
    V_ckT = 0.*(kf_ck.*(PCrT).*ADPT-kr_ck.*CrT.*ATPT);
    
    ATP_inhT = power((1.00000+nATP.*(ATPT./KI_ATP))./(1.00000+ATPT./KI_ATP), 4.00000);
    V_pkT = enfatraT.*Vmax_pk.*(GAPT./(GAPT+Km_GAP_pk)).*(ADPT./(ADPT+Km_ADP_pk)).*ATP_inhT;
    pa=1;
    V_opT = enfamitT.*Vmax_op.*((pa.*PYRT)./((pa.*PYRT)+Km_PYR_op)).*(ADPT./(ADPT+Km_ADP_op)).*(1.00000./(1.00000+0.100000.*(ATPT./ADPT)));
    
    AMPT = ANP-(ATPT+ADPT);
    AMP_actT = power((1.00000+AMPT./Ka_AMP)./(1.00000+nAMP.*(AMPT./Ka_AMP)), 4.00000);
    V_pfkT = Vmax_pfk.*(F6PT./(F6PT+Km_F6P_pfk)).*(ATPT./(ATPT+Km_ATP_pfk)).*(F26PT./(F26PT+Km_F26P_pfk)).*ATP_inhT.*AMP_actT;
    
    eta_ldhT=1-beta_ldh_ros.*(1./(1+power(Kldh_ros./ROST,4)));%((ROST.^4)./((ROST.^4)+(Kldh_ros.^4)));
    V_ldhT = 1.*eta_ldhT.*(kf_ldh.*PYRT-kr_ldh.*LACT);
    V_lacT = (Vlac_0.*(1.00000+v_stim1.*K_lac)-K_lac_eff.*LACT);
    
    V_hkT = Vmax_hk.*(ATPT./(ATPT+Km_ATP_hk)).*(power(1.00000+power(F6PT./KI_F6P, 4.00000), -1.00000)).*GLCe;
    AMP_pfk2T = (power(AMPT./Kamp_pfk2, nh_amp))./(1.00000+power(AMPT./Kamp_pfk2, nh_amp));
    V_pfk2T = Vmaxf_pfk2.*(ATPT./(ATPT+Km_ATP_pfk2)).*(F6PT./(F6PT+Km_F6P_pfk2)).*AMP_pfk2T-Vmaxr_pfk2.*(F26PT./(F26PT+Km_F26P_pfk2));
    
    V_ATPaseT = Vmax_ATPase.*(ATPT./(ATPT+Km_ATP)).*(1.00000+v_stim);
    dAMP_dATPT = -1.00000+Q_adk./2.00000+-(0.500000.*(power(uADPT, 1/2)))+Q_adk.*(ANP./(ATPT.*(power(uADPT, 1/2))));
    
    GSSGT=(GSH_totT-GSHT)./2;
    NADPT=NADPH_tot-NADPHT;
    VpppT = Vmax_ppp.*(F6PT./(F6PT+Km_F6P_pfk))./(1+((NADPHT./NADPT)./Ki_nadph));
    VgrT = kf_gr.*GSSGT.*NADPHT-kr_gr.*GSHT.*NADPT;
    
    % PD pathology pathways
    eta_opT=eta_op_max-beta_op_asyn.*(1./(1+power(Kasyn./ASYNAT,4)));%(((ASYNAT.^4)./((ASYNAT.^4)+(Kasyn.^4))));
    
    Vros_leakT=(0.5282./ATPT).*(1-eta_opT).*V_opT;%0.221 or (0.5282./ATP)
    %     Vros_leak=0.221.*(1-eta_op).*V_op;%0.221 or (0.5282./ATP)
    
    Vros_catT=Kros_cat.*ROST;
    
    Vros_dopaT=0.001.*Kros_dopa.*(((cdaT)./((cdaT)+(Kasyn.*1000))));
    
    Vros_doxT=Kros_dox.*GSHT.*ROST;
    
    Vasyn_synT=Kasyn_syn;
    
    Vasyn_oxT=Kasyn_ox.*ROST.*ASYN_T;
    
    Vasyn_toT=Kasyn_to.*ASYN_T;
    
    Vasyn_aggT=Krasyn_agg.*ASYNAT.*(1./(1+power(Kasyn_agg./ASYNAT,6)));%(((ASYNAT.^6)./((ASYNAT.^6)+(Kasyn_agg.^6))));
    
    UbT=Ub_tot-ASYNTT;
    Vasyn_tagT=Kasyn_tag.*ASYNAT.*UbT.*ATPT;
    
    Vasyn_prtT=Krasyn_prt.*ASYNTT.*ATPT.*(1-beta_asyn_prt.*(1./(1+power(Kasyn_prt./ASYNGT,4))));%(((ASYNGT.^4)./((ASYNGT.^4)+(Kasyn_prt.^4)))));
    
    Vasyn_lysoT=Kasyn_lyso.*ASYNGT.*ATPT;
    
    Vasyn_lbT=Krasyn_lb.*ASYNGT.*(1./(1+power(Kasyn_lb./ASYNGT,6)));%(((ASYNGT.^6)./((ASYNGT.^6)+(Kasyn_lb.^6))));
    
    extT=(power(1-dAMP_dATPT, -1.00000));
    extT((isnan(extT)| isinf(extT)))=0;
    
    vAADCextT = 0.24.*(vAADC_Vmax)./(1+(power((vAADC_Km./sLDT),0.5)));
    
    %%
    %%%%%%%%%%%%%%%%%%%%%% Differential equations %%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    cdanxtT = cdaT+(0*jsyntT + jdatT - jvmatT - jidaT + jldopaT).*dt;%cda
    vdanxtT = vdaT+(jvmatT - jrelT).*dt;%vda
    edanxtT = edaT+(jrelT - jdatT - jedaT  + vAADCextT).*dt;%eda
    F6PnxtT = F6PT+(V_hkT-(V_pfkT-V_pfk2T)-VpppT.*(1./6)).*dt;%F6P
    F26PnxtT = F26PT+(V_pfk2T).*dt;%F26P
    GAPnxtT = GAPT+(2.00000.*V_pfkT-V_pkT).*dt;%GAP
    PYRnxtT = PYRT+(V_pkT-(V_opT+V_ldhT)).*dt;%PYR
    LACnxtT = LACT+(2.25000.*V_ldhT+V_lacT).*dt;%LAC
    ATPnxtT = ATPT+(((1.*(2.*V_pkT+15.*eta_opT.*V_opT+V_ckT))-(V_hkT+V_pfkT+V_pfk2T+V_ATPaseT+V_pumpsT+25.*Vasyn_prtT+1.*(3.*Vasyn_tagT+10.*Vasyn_lysoT))).*extT).*dt;%ATP
    PCrnxtT = PCrT+(-V_ckT).*dt;%PCr
    ROSnxtT = ROST+(Vros_leakT + Vros_exT - Vros_catT + Vros_dopaT - Vros_doxT).*dt;%ROS
    ASYNnxtT = ASYN_T+(Vasyn_synT - Vasyn_oxT - Vasyn_toT).*dt;%ASYN
    ASYNAnxtT = ASYNAT+(Vasyn_oxT - Vasyn_aggT - Vasyn_tagT).*dt;%ASYNA
    ASYNTnxtT = ASYNTT+(Vasyn_tagT - Vasyn_prtT).*dt;%ASYNT
    ASYNGnxtT = ASYNGT+(Vasyn_aggT - Vasyn_lysoT - Vasyn_lbT).*dt;%ASYNG
    LBnxtT = LBT+(Vasyn_lbT).*dt;%LB
    NADPHnxtT = NADPHT+(2.*VpppT - VgrT).*dt;%NADPH
    GSHnxtT = GSHT+(2.*VgrT - 2.*Vros_doxT).*dt;%
    LDOPAnxtT = LDOPAT+(((Vtran_max.*sLDT)./(Ksld.*(1+(sTYR./Kstyr)+(sTRP./Kstrp))+sLDT)) - jldopaT + jsyntT).*dt;%ldopa
    
    cdaT=cdanxtT;vdaT=vdanxtT;edaT=edanxtT;
    NADPHT=NADPHnxtT;
    GSHT=GSHnxtT;
    F6PT=F6PnxtT;
    F26PT=F26PnxtT;
    GAPT=GAPnxtT;
    PYRT=PYRnxtT;
    LACT=LACnxtT;
    ATPT=ATPnxtT;
    PCrT=PCrnxtT;
    ROST=ROSnxtT;
    ASYN_T=ASYNnxtT;
    ASYNAT=ASYNAnxtT;
    ASYNTT=ASYNTnxtT;
    ASYNGT=ASYNGnxtT;
    LBT=LBnxtT;
    LDOPAT=LDOPAnxtT;
    
    ATPT(isnan(ATPT)| isinf(ATPT)| ATPT<0)=0;
    ATP(isnan(ATP)| isinf(ATP)| ATP<0)=0;
    cdaT(isnan(cdaT)| isinf(cdaT)| cdaT<0)=0;
    eda(isnan(eda)| isinf(eda)| eda<0)=0;
    edaT(isnan(edaT)| isinf(edaT)| edaT<0)=0;
    ROST(isnan(ROST)| isinf(ROST)| ROST<0)=0;
    
    % Killing based on ROS
    indsapT=find(ROST>rosthr);
    indsappT=[indsappT indsapT'];
    if isempty(indsapT)==0
        indsappT=unique(indsappT,'stable');
        indsapT=[];
    end
    
    Ca_iT(indsappT)=1.*zeros(size(indsappT));
    cdaT(indsappT) = 1e-4.*zeros(size(indsappT)); %mM
    vdaT(indsappT) = 500.*zeros(size(indsappT)); %mM
    edaT(indsappT) = 26e-6.*zeros(size(indsappT)); %mM
    ATPT(indsappT)=0.*zeros(size(indsappT));
    
    %     inds=find(V_snc_array(k) >= -20 && V_snc_array(k) > V_snc_array(k-1) && V_snc_array(k) > V_snc_array(k+1));
    %     inds=find(V_snc <=80 & V_snc >=-20);
    %     snc_firings=[snc_firings; k+0*inds,inds+0*inds];
    
    %----------------------------------------STN-----------------------------------------%
    %---------------------------Input from stn to stn(laterals)--------------------------%
    
    % psp variable
    h_nmdastn = (1-lam_nmda).* h_nmdastn + lam_nmda.*stncurr_spk; %psp nmda stn lat
    h_ampastn = (1-lam_ampa).* h_ampastn + lam_ampa.*stncurr_spk; %psp ampa stn lat
    
    tmplat_nmda_stn = h_nmdastn.*(Enmda - Vstn);
    tmplat_ampa_stn = h_ampastn.*(Eampa - Vstn);
    
    Ilat_nmda_stn = conv2(tmplat_nmda_stn, wlatstn, 'same'); % lat nmda stn
    Ilat_ampa_stn = conv2(tmplat_ampa_stn, wlatstn, 'same');% lat ampa stn
    
    B = 1./(1 + (mg0./3.57).*exp(-0.062.*Vstn));
    
    %----------------------------------------Input from gpe to stn---------------------%
    % psp variable
    h_gs = (1-lam_gaba).* h_gs + lam_gaba.*gpecurr_spk;% input from gpe to nmda stn
    % gaba current
    I_gs = wgs.*h_gs.*(Egaba - Vstn);
    
    Istnstn=B.*Ilat_nmda_stn + Ilat_ampa_stn;
    % total current stn recieves
    Itmpstn =  Istnstn + I_gs; % total currents
    
    % V,U updated
    dvstn = taustn.*(((0.04.*Vstn.*Vstn)+5.*Vstn - Ustn +140 + Istn + Itmpstn)./Cstn);%+wgn(Mstn,Nstn,lnoise_st,imp_st);
    dustn = taustn.*astn.*(bstn.*(Vstn) - Ustn);
    Vstn_nxt = Vstn + dvstn;
    Ustn_nxt = Ustn + dustn;
    
    inds = find(Vstn_nxt > vpeak_stn);
    %     stn_firings=[stn_firings; k+0*inds,inds];
    
    Vstn_nxt(inds) = cstn.*ones(size(inds));
    Ustn_nxt(inds) = Ustn(inds) + dstn.*ones(size(inds));
    stncurr_spk = stn_zeros;
    stncurr_spk(inds) = ones(size(inds));
    
    Vstn = Vstn_nxt;
    Ustn = Ustn_nxt;
    
    %----------------------------------------GPe-----------------------------------------%
    %--------------------- Input from stn to gpe----------------------------------%
    % psp variable
    h_nmdagpe = (1-lam_nmda).* h_nmdagpe + lam_nmda.*stncurr_spk;
    h_ampagpe = (1-lam_ampa).* h_ampagpe + lam_ampa.*stncurr_spk;
    
    % nmda and ampa currents
    Inmdagpe = wsg.*h_nmdagpe.*(Enmda - Vgpe);
    Iampagpe = wsg.*h_ampagpe.*(Eampa - Vgpe);
    
    %-------------------------- Input from gpe to gpe(laterals)-------------------%
    % psp variable
    hlat_gaba_gpe = (1-lam_gaba).* hlat_gaba_gpe + lam_gaba.*gpecurr_spk;
    tmplat_gaba_gpe = wda_gpe.*hlat_gaba_gpe.*(Egaba - Vgpe);
    
    % lateral currents
    Ilatgpe = conv2(tmplat_gaba_gpe, wlatgpe, 'same');
    
    B = 1./(1 + (mg0/3.57).*exp(-0.062.*Vgpe));
    I_sg=B.*Inmdagpe+Iampagpe;
    Itmpgpe = Ilatgpe + I_sg;
    
    % LFP
    % LFP_GPe_exc(k)=sum(sum(I_sg));
    % LFP_GPe_inh(k)=sum(sum(Ilatgpe));
    % LFP_GPe_tot(k)=LFP_GPe_exc(k)+LFP_GPe_inh(k);
    
    % V ,U updated
    dvgpe = taugpe.*((0.04.*Vgpe.*Vgpe)+5.*Vgpe+140 - Ugpe + Itmpgpe+Igpe)./Cgpe;
    dugpe = taugpe.*agpe.*(bgpe.*(Vgpe) - Ugpe);
    Vgpe_nxt = Vgpe + dvgpe;
    Ugpe_nxt = Ugpe + dugpe;
    
    inds = find(Vgpe_nxt > vpeak_gpe);
    %     gpe_firings=[gpe_firings; k+0*inds,inds];
    
    Vgpe_nxt(inds) = cgpe.*ones(size(inds));
    Ugpe_nxt(inds) = Ugpe(inds) + dgpe.*ones(size(inds));
    gpecurr_spk = gpe_zeros;
    gpecurr_spk(inds) = ones(size(inds));
    
    Vgpe = Vgpe_nxt;
    Ugpe = Ugpe_nxt;
    
    %----------------------------------------LIT-----------------------------------------%
    %     terDA=zeros(1,PsncT);
    %     for nn=1:Pmsn
    %         terDA(nn)=mean(edaT(idx_DA_MSN(nn,:)));
    %     end
    
    terDA=mean(edaT(idx_DA_MSN(:,:)),2);
    terDA=reshape(terDA,Mmsn,Nmsn);
    
    %     phi1=noDAT.*RescaleRange(terDA,1e-6,2e-6,0,1);%1-3
    %     phi1=RescaleRange(mean2(edaT),1e-6,5e-6,0,1);
    %     phi1=1.*RescaleRange(terDA/(0.005*k),1e-6,3e-5,0,1);
    
    if pdon==1
        phi1=noDAT.*RescaleRange(terDA/(0.005*k),1e-6,2e-6,0,1);
    else
        phi1=noDAT.*RescaleRange(terDA,1e-6,2e-6,0,1);
    end
    
    phi1(phi1>1)=1;
    phi1(phi1<0)=0;
    
    phi1(isnan(phi1)| isinf(phi1)| phi1<0)=0;
    
    %     if phi1<0
    %         phi1=0;
    %     end
    %     if phi1>1
    %         phi1=1;
    %     end
    
    
    
    %---------------------------------------CTR---------------------------------------%
    
    % V,U updated
    dvctr = tauctr.*((0.7.*(Vctr+60).*(Vctr+40) - Uctr + Ictr)./Cctr);%+wgn(Mctr,Nctr,lnoise_st,imp_st);
    ductr = tauctr.*actr.*(bctr.*(Vctr+60) - Uctr);
    Vctr_nxt = Vctr + dvctr;
    Uctr_nxt = Uctr + ductr;
    
    inds = find(Vctr_nxt > vpeak_ctr);
    %     ctr_firings=[ctr_firings; k+0*inds,inds];
    
    Vctr_nxt(inds) = cctr.*ones(size(inds));
    Uctr_nxt(inds) = Uctr(inds) + dctr.*ones(size(inds));
    ctrcurr_spk = ctr_zeros;
    ctrcurr_spk(inds) = ones(size(inds));
    
    Vctr = Vctr_nxt;
    Uctr = Uctr_nxt;
    
    %---------------------------------------MSN---------------------------------------%
    
    Vmsn1r=Vmsn11r.*(1+(Kmsn1.*phi1));
    dmsn1=dmsn11.*(1-(Lmsn1.*phi1));
    %     Iextmsn1=350*(1+beta1*phi1);
    %     Imsn1=Iextmsn1;
    
    Vmsn2r=Vmsn22r.*(1+(Kmsn2.*(phi1./2)));
    dmsn2=dmsn22.*(1-(Lmsn2.*(phi1./2)));
    
    h_nmda_msn1 = (1-lam_nmda).* h_nmda_msn1 + lam_nmda.*ctrcurr_spk;
    h_ampa_msn1 = (1-lam_ampa).* h_ampa_msn1 + lam_ampa.*ctrcurr_spk;
    
    % nmda and ampa currents
    Inmda_msn1 = wctr_msn1.*h_nmda_msn1.*(Enmda - Vmsn1);
    Iampa_msn1 = wctr_msn1.*h_ampa_msn1.*(Eampa - Vmsn1);
    B11 = 1./(1 + (mg0./3.57).*exp(-0.062.*Vmsn1));
    Imsn1 = (Iampa_msn1 + (B11.*Inmda_msn1))*(1+beta1.*phi1);
    %     Ictr_msn1=[Ictr_msn1 Imsn1];
    
    % psp variable for D12
    h_nmda_msn2 = (1-lam_nmda).* h_nmda_msn2 + lam_nmda.*ctrcurr_spk;
    h_ampa_msn2 = (1-lam_ampa).* h_ampa_msn2 + lam_ampa.*ctrcurr_spk;
    
    % nmda and ampa currents
    Inmda_msn2 = wctr_msn2.*h_nmda_msn2.*(Enmda - Vmsn2);
    Iampa_msn2 = wctr_msn2.*h_ampa_msn2.*(Eampa - Vmsn2);
    B12 = 1./(1 + (mg0./3.57).*exp(-0.062.*Vmsn2));
    Imsn2 = (Iampa_msn2 + (B12.*Inmda_msn2))*(1+beta1.*(phi1./2));
    %     Ictr_msn2=[Ictr_msn2 Imsn2];
    
    %----------------------------------------msn1-----------------------------------------%
    
    % MSN1 to MSN2------%
    hlat_gaba_msn = (1-lam_gaba).* hlat_gaba_msn + lam_gaba.*msn1curr_spk;
    I_gaba_msn = w_msn1_msn2*(hlat_gaba_msn.*(Egaba - Vmsn2));
    %     I_msn1_msn2=[I_msn1_msn2 I_gaba_msn];
    
    % V,U updated
    dvmsn1 = taumsn1.*((kmsn1.*((Vmsn1-Vmsn1r).*(Vmsn1-Vmsn1t))- Umsn1 + Imsn1)./Cmsn1);%+wgn(Mmsn1,Nmsn1,lnoise,imp);
    dumsn1 = taumsn1.*amsn1.*(bmsn1.*(Vmsn1-Vmsn1r) - Umsn1);
    %     dvmsn1 = taumsn1.*(((0.04.*Vmsn1.*Vmsn1)+5.*Vmsn1 - Umsn1 + Itmpmsn1 +140)./Cmsn1)+sqrt(taumsn1)*(lnoise*randn);
    %     dumsn1 = taumsn1.*amsn1.*(bmsn1.*(Vmsn1) - Umsn1);
    Vmsn1_nxt = Vmsn1 + dvmsn1;
    Umsn1_nxt = Umsn1 + dumsn1;
    
    inds = find(Vmsn1_nxt > vpeak_msn1);
    %     msn1_firings=[msn1_firings; k+0*inds,inds];
    
    Vmsn1_nxt(inds) = cmsn1.*ones(size(inds));
    Umsn1_nxt(inds) = Umsn1(inds) + dmsn1(inds).*ones(size(inds));
    msn1curr_spk = msn1_zeros;
    msn1curr_spk(inds) = ones(size(inds));
    
    
    Vmsn1 = Vmsn1_nxt;
    Umsn1 = Umsn1_nxt;
    
    %     dVmsn11(k)=Vmsn1;
    % dVsnc1(k)=V_snc;
    
    %----------------------------------------msn2-----------------------------------------%
    
    % V,U updated
    dvmsn2 = taumsn2.*((kmsn2.*((Vmsn2-Vmsn2r).*(Vmsn2-Vmsn2t))- Umsn2 +Imsn2+1.*I_gaba_msn)./Cmsn2);%+wgn(Mmsn2,Nmsn2,lnoise,imp);
    dumsn2 = taumsn2.*amsn2.*(bmsn2.*(Vmsn2-Vmsn2r) - Umsn2);
    %     dvmsn2 = taumsn2.*(((0.04.*Vmsn2.*Vmsn2)+5.*Vmsn2 - Umsn2 + Itmpmsn2 +140)./Cmsn2)+sqrt(taumsn2)*(lnoise*randn);
    %     dumsn2 = taumsn2.*amsn2.*(bmsn2.*(Vmsn2) - Umsn2);
    Vmsn2_nxt = Vmsn2 + dvmsn2;
    Umsn2_nxt = Umsn2 + dumsn2;
    
    inds = find(Vmsn2_nxt > vpeak_msn2);
    %     msn2_firings=[msn2_firings; k+0*inds,inds];
    
    Vmsn2_nxt(inds) = cmsn2.*ones(size(inds));
    Umsn2_nxt(inds) = Umsn2(inds) + dmsn2(inds).*ones(size(inds));
    msn2curr_spk = msn2_zeros;
    msn2curr_spk(inds) = ones(size(inds));
    
    Vmsn2 = Vmsn2_nxt;
    Umsn2 = Umsn2_nxt;
    
    % Killing based on apoptosis threshold
    indsap=find(apop>apopthr);
    indsapp=[indsapp indsap'];
    if isempty(indsap)==0
        indsapp=unique(indsapp,'stable');
        indsap=[];
    end
    V_snc(indsapp) = -80.*ones(size(indsapp));
    apop(indsapp)=zeros(size(indsapp));
    Ca_i(indsapp)=zeros(size(indsapp));
    cda(indsapp) = 1e-4.*zeros(size(indsapp)); %mM
    vda(indsapp) = 500.*zeros(size(indsapp)); %mM
    eda(indsapp) = 26e-6.*zeros(size(indsapp)); %mM
    ATPused(indsapp)=0.*zeros(size(indsapp));
    ATP(indsapp)=0.*zeros(size(indsapp));
    
    % Storing pattern of SNc cell loss
    %     indsnum = find(V_snc < -70);
    %     if isempty(indsnum)==0
    %         indsnums=setdiff(indsnum,indsappp);
    %         if isempty(indsnums)==0
    %             indsappp=[indsappp; k+0*indsnums,indsnums];
    %         end
    %     end
    %
    
    %     clit{k}=setdiff(indsapp,indsappp);
    
    %     indsappp=indsapp;
    
    eda1=sum(sum(eda))/(Psnc);
    
    yeda1(count)=sum(sum(eda))/(Psnc);
    yatp(count)=sum(sum(ATP))/(Psnc);
    
    %     dDA1(count)=DA;
    dnc1(count)=numel(indsapp);
    dnc2(count)=numel(indsappT);
    
    Inc=(numel(indsapp));
    IncT=(numel(indsappT));
    
    yeda1T(count)=sum(sum(edaT))./(PsncT);
    ycda1T(count)=sum(sum(cdaT))./(PsncT);
    yatpT(count)=sum(sum(ATPT))./(PsncT);
    yDAT(count)=sum(sum(phi1))./PsncT;
    Isp(count)=sum(sum(Isp_snc))./Psnc;
    
    ycaiS(count)=sum(sum(Ca_i))./(Psnc);
    ycaiT(count)=sum(sum(Ca_iT))./(PsncT);
    
    ycamtS(count)=sum(sum(Ca_mt))./(PsncT);
    yrosT(count)=sum(sum(ROST))./(PsncT);
    %     yrosT(:,count)=reshape(ROST,PsncT,1);
    %     ycamtS(:,count)=reshape(Ca_mt,Psnc,1);
    
    %     if yatp(count) <= 0
    %         bre=1;
    %         break
    %     end
    
    count=count+1;
    % Sub-sampling iterative
    if (k==samp_start)
        
        if gpuon==1
            yeda0=sub_sampling_GPU(yeda1,dt);
            yatp0=sub_sampling_GPU(yatp,dt);
            yeda0T=sub_sampling_GPU(yeda1T,dt);
            ycda0T=sub_sampling_GPU(ycda1T,dt);
            yatp0T=sub_sampling_GPU(yatpT,dt);
            yros0T=sub_sampling_GPU(yrosT,dt);
            dnc0=sub_sampling_GPU(dnc1,dt);
            dnc02=sub_sampling_GPU(dnc2,dt);
            yDA0T=sub_sampling_GPU(yDAT,dt);
            Isp0=sub_sampling_GPU(Isp,dt);
            ycai0S=sub_sampling_GPU(ycaiS,dt);
            ycai0T=sub_sampling_GPU(ycaiT,dt);
            ycamt0S=sub_sampling_GPU(ycamtS,dt);
            %             dDA0=sub_sampling_GPU(dDA1,dt);
        elseif gpuon==0
            yeda0=sub_sampling(yeda1,dt);
            yatp0=sub_sampling(yatp,dt);
            yeda0T=sub_sampling(yeda1T,dt);
            ycda0T=sub_sampling(ycda1T,dt);
            yatp0T=sub_sampling(yatpT,dt);
            yros0T=sub_sampling(yrosT,dt);
            dnc0=sub_sampling(dnc1,dt);
            dnc02=sub_sampling(dnc2,dt);
            yDA0T=sub_sampling(yDAT,dt);
            Isp0=sub_sampling(Isp,dt);
            ycai0S=sub_sampling(ycaiS,dt);
            ycai0T=sub_sampling(ycaiT,dt);
            ycamt0S=sub_sampling(ycamtS,dt);
            %             dDA0=sub_sampling(dDA1,dt);
        end
        
        deda=[deda yeda0];
        datp=[datp yatp0];
        dedaT=[dedaT yeda0T];
        dcdaT=[dcdaT ycda0T];
        datpT=[datpT yatp0T];
        drosT=[drosT yros0T];
        dncS=[dncS dnc0];
        dncT=[dncT dnc02];
        dDAT=[dDAT yDA0T];
        dIsp=[dIsp Isp0];
        dcaiS=[dcaiS ycai0S];
        dcaiT=[dcaiT ycai0T];
        dcamtS=[dcamtS ycamt0S];
        %         dDA=[dDA dDA0];
        
        samp_start=samp_start+(stt/dt);
        count=1;
        disp(k*dt)
    end
    %     dVstn=[dVstn Vstn(16,16)];
    %     dVgpe=[dVgpe Vgpe(16,16)];
    %     dVsnc=[dVsnc V_snc(4,4)];
    %     dVctr=[dVctr Vctr(16,16)];
    %     dVmsn1=[dVmsn1 Vmsn1(4,4)];
    %     dVmsn2=[dVmsn2 Vmsn2(16,16)];
    
    
    %     disp(k*dt)
    
end

if bre==1
    f10=strcat('Break_',filename);
    save(f10);
end

% kid=((Psnc)-dnc);
% out_cell={};
% out_cell{1,1}=clit;
% out_cell{2,1}=kid;
% out_cell{3,1}=deda;
% out_cell{4,1}=dDA;

% [snc_firings]=ConvertAPtoST(snc_firings,Psnc);

% base0=Pmsn/2;
% msn1frequency=size(msn1_firings,1)/(2*base0.*dt.*(numel(tspan)).*1e-3);
%
% base0=Pmsn/2;
% msn2frequency=size(msn2_firings,1)/(2*base0.*dt.*(numel(tspan)).*1e-3);
%
% base0=Psnc/2;
% sncfrequency=size(snc_firings,1)/(2*base0.*dt.*(numel(tspan)).*1e-3);
%
% base0=Pstn/2;
% stnfrequency=size(stn_firings,1)/(2*base0.*dt.*(numel(tspan)).*1e-3);
%
% base0=Pgpe/2;
% gpefrequency=size(gpe_firings,1)/(2*base0.*dt.*(numel(tspan)).*1e-3);
%
% base0=Pctr/2;
% ctrfrequency=size(ctr_firings,1)/(2*base0.*dt.*(numel(tspan)).*1e-3);

stt=toc;
if stt>60 && stt<=3600
    stt1=stt/60;
    f1=['Simulation time is ',num2str(stt1),' minutes_',num2str(sLD)];
    %     disp(f1);
elseif stt>3600
    stt2=stt/(60*60);
    f1=['Simulation time is ',num2str(stt2),' hours_',num2str(sLD)];
    %     disp(f1);
else
    f1=['Simulation time is ',num2str(stt),' seconds_',num2str(sLD)];
    %     disp(f1);
end
simtime=f1;