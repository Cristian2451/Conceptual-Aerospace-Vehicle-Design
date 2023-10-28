clear
clc

W0S = linspace(0,100000,10000);

Clmax = 2;
Clmax_1 = 2.45;
Clmax_2 = 3;

%%%%%%%%%%%%%%%%%%%%       Cruise       %%%%%%%%%%%%%%%%%%%%%

alpha_cruise_1 = 0.948;
alpha_cruise_2 = 0.874;
alpha_cruise_3 = 0.814;
alpha_cruise_missed = 0.783;
alpha_loit = 0.765; 

rho_1 = 0.5489; %25000ft
rho_2 = 0.52045; %26500ft
rho_3 = 0.4842; %28500ft
rho_missed = 0.4754; %29000ft
rho_loit = 1.056; %5000ft

% beta_cruise_1 = 0.4; %25000ft 0.69M
% beta_cruise_2 = 0.35; %29000ft 0.75M

beta_cruise_1 = 0.29; %25000ft 0.69M
beta_cruise_2 = 0.25; %29500ft 0.75M
beta_cruise_3 = 0.24;
beta_cruise_4 = 0.23;
beta_loit = 0.525;
%beta_loit = 0.7;

M_1 = 0.69;
M_2 = 0.75;

V_tmd = 117.042; % V_imd in m/s for loiter
V_inf_1 = M_1*309.7; %TAS Cruise 1
V_inf_2 = M_2*303.5;
V_inf_3 = M_2*300.4;
V_inf_missed = M_2*299;

e = 1/1.2;
AR = 10;
n = 1;


LDcruise = 16.44;
LDloiter = 20;
LDmax = 20;
Cd0 = pi*AR*e/(2*LDmax)^2;

% Cd0_cruise = pi*AR*e/(2*LDcruise)^2;
% Cd0_loit = pi*AR*e/(2*LDloiter)^2;

TW0_cruise_1 = (alpha_cruise_1/beta_cruise_1)*((0.5*rho_1*V_inf_1^2*Cd0./(alpha_cruise_1.*W0S))+(alpha_cruise_1*n^2.*W0S/(0.5*rho_1*V_inf_1^2*pi*AR*e)));
TW0_cruise_2 = (alpha_cruise_2/beta_cruise_2)*((0.5*rho_2*V_inf_2^2*Cd0./(alpha_cruise_2.*W0S))+(alpha_cruise_2*n^2.*W0S/(0.5*rho_2*V_inf_2^2*pi*AR*e)));
TW0_cruise_3 = (alpha_cruise_3/beta_cruise_3)*((0.5*rho_3*V_inf_3^2*Cd0./(alpha_cruise_3.*W0S))+(alpha_cruise_3*n^2.*W0S/(0.5*rho_3*V_inf_3^2*pi*AR*e)));
TW0_cruise_missed = (alpha_cruise_missed/beta_cruise_4)*((0.5*rho_missed*V_inf_missed^2*Cd0./(alpha_cruise_missed.*W0S))+(alpha_cruise_missed*n^2.*W0S/(0.5*rho_missed*V_inf_missed^2*pi*AR*e)));


%%%%%%%%%%%%%%%%%%%%       Take-off       %%%%%%%%%%%%%%%%%%%%%

sigma = 1;
Ne = 2;
beta = 1;
alpha_TO_1 = 0.98;
alpha_TO_2 = 0.915;
alpha_TO_3 = 0.843;

TODA_1 = 2300*(100/115);
TODA_2 = 3288*(100/115);
TODA_3 = 2119*(100/115);

%ClLOF = 0.8*Clmax;

TW0_TO_1 = (alpha_TO_1/beta)*(0.297-0.019*Ne).*W0S/(sigma*0.8*Clmax*TODA_1);
TW0_TO_2 = (alpha_TO_2/beta)*(0.297-0.019*Ne).*W0S/(sigma*0.8*Clmax*TODA_2);
TW0_TO_3 = (alpha_TO_3/beta)*(0.297-0.019*Ne).*W0S/(sigma*0.8*Clmax*TODA_3);

TW0_TO_11 = (alpha_TO_1/beta)*(0.297-0.019*Ne).*W0S/(sigma*0.8*Clmax_1*TODA_1);
TW0_TO_12 = (alpha_TO_2/beta)*(0.297-0.019*Ne).*W0S/(sigma*0.8*Clmax_1*TODA_2);
TW0_TO_13 = (alpha_TO_3/beta)*(0.297-0.019*Ne).*W0S/(sigma*0.8*Clmax_1*TODA_3);

TW0_TO_21 = (alpha_TO_1/beta)*(0.297-0.019*Ne).*W0S/(sigma*0.8*Clmax_2*TODA_1);
TW0_TO_22 = (alpha_TO_2/beta)*(0.297-0.019*Ne).*W0S/(sigma*0.8*Clmax_2*TODA_2);
TW0_TO_23 = (alpha_TO_3/beta)*(0.297-0.019*Ne).*W0S/(sigma*0.8*Clmax_2*TODA_3);


%%%%%%%%%%%%%%%%%%%%       Landing       %%%%%%%%%%%%%%%%%%%%%

Sa = 305;
%Kr = 0.66;
sigma = 1;

LDA_1 = 1815*0.6;
LDA_2 = 3288*0.6;
LDA_3 = 2059*0.6;

W0S_emergency = (1/0.975)*(LDA_1 - Sa)*sigma*Clmax/(0.51*0.77);
W0S_1 = (1/0.939)*(LDA_2 - Sa)*sigma*Clmax/(0.51*0.77);
W0S_2 = (1/0.866)*(LDA_3 - Sa)*sigma*Clmax/(0.51*0.77);
W0S_3 = (1/0.765)*(LDA_1 - Sa)*sigma*Clmax/(0.51*0.77);

W0S_emergency_1 = (1/0.975)*(LDA_1 - Sa)*sigma*Clmax_1/(0.51*0.77);
W0S_11 = (1/0.939)*(LDA_2 - Sa)*sigma*Clmax_1/(0.51*0.77);
W0S_12 = (1/0.866)*(LDA_3 - Sa)*sigma*Clmax_1/(0.51*0.77);
W0S_13 = (1/0.765)*(LDA_1 - Sa)*sigma*Clmax_1/(0.51*0.77);

W0S_emergency_2 = (1/0.975)*(LDA_1 - Sa)*sigma*Clmax_2/(0.51*0.77);
W0S_21 = (1/0.939)*(LDA_2 - Sa)*sigma*Clmax_2/(0.51*0.77);
W0S_22 = (1/0.866)*(LDA_3 - Sa)*sigma*Clmax_2/(0.51*0.77);
W0S_23 = (1/0.765)*(LDA_1 - Sa)*sigma*Clmax_2/(0.51*0.77);

%%%%%%%%%%%%%%%%%%%%       Stall       %%%%%%%%%%%%%%%%%%%%%

% CLmax=1.5;
% 
% V_stall=V_inf_2/1.2;
% W0S_stall=0.5*rho_2*V_stall^2*CLmax;

%%%%%%%%%%%%%%%%%%%%       Climb       %%%%%%%%%%%%%%%%%%%%%

% dhdt_climb_1=12.7;
% climb_grad=2.4; % minimum
% climb_angle=atand(climb_grad/100);
% TW0_climb=(Ne/(Ne-1))*(1/15 +sind(climb_angle));

% FAR climb requirements for OEI
CGR_FSC = sind(atand(0/100));
LD_FSC = 0.5*sqrt(pi*(e-0.05)*AR/(Cd0+0.02+0.02));
first_segment_climb = (Ne/(Ne-1))*(1/(LD_FSC) + CGR_FSC);

CGR_SSC = sind(atand(2.4/100));
LD_SSC = 0.5*sqrt(pi*(e-0.05)*AR/(Cd0+0.02));
second_segment_climb = (Ne/(Ne-1))*(1/(LD_SSC) + CGR_SSC);

CGR_TSC = sind(atand(1.2/100));
LD_TSC = 0.5*sqrt(pi*(e)*AR/(Cd0));
third_segment_climb = (Ne/(Ne-1))*(1/(LD_TSC) + CGR_TSC);

CGR_AG = sind(atand(2.1/100));
LD_AG = 0.5*sqrt(pi*(e-0.06)*AR/(Cd0+0.045));
approach_goaround = (Ne/(Ne-1))*(1/(LD_AG) + CGR_AG);

CGR_LG = sind(atand(3.2/100));
LD_LG = 0.5*sqrt(pi*(e-0.1)*AR/(Cd0+0.02+0.07));
landing_goaround = (1/(LD_LG) + CGR_LG);


%%%%%%%%%%%%%%%%%%%%       Cruise Turn      %%%%%%%%%%%%%%%%%%%%%

TW0_cruise_turn = (alpha_cruise_1/beta_cruise_1)*(0.5*rho_1*V_inf_1^2*Cd0./(alpha_cruise_1.*W0S)+alpha_cruise_1*(1/cosd(12.5))^2.*W0S/(0.5*rho_1*V_inf_1^2*3.141*AR*e));
TW0_loiter_turn = (alpha_cruise_1/beta_loit)*(0.5*rho_loit*V_tmd^2*Cd0./(alpha_cruise_1.*W0S)+alpha_cruise_1*(1/cosd(25))^2.*W0S/(0.5*rho_loit*V_tmd^2*3.141*AR*e));

%%%%%%%%%%%%%%%%%%%%       Cruise M=0.75 (spec point)       %%%%%%%%%%%%%%%%%%%%%

Mmax = 0.75;
V_inf_1max = Mmax*309.7;


TW0_cruise_1max = (alpha_cruise_1/beta_cruise_1)*((0.5*rho_1*V_inf_1max^2*Cd0./(alpha_cruise_1.*W0S))+(alpha_cruise_1*n^2.*W0S/(0.5*rho_1*V_inf_1max^2*pi*AR*e)));
%TW0_cruise_2max = (alpha_cruise_2/beta_cruise_2)*((0.5*rho_2*V_inf_2max^2*Cd0./(alpha_cruise_2.*W0S))+(alpha_cruise_2*n^2.*W0S/(0.5*rho_2*V_inf_2max^2*pi*AR*e)));
%TW0_cruise_3max = (alpha_cruise_3/beta_cruise_2)*((0.5*rho_3*V_inf_3max^2*Cd0./(alpha_cruise_3.*W0S))+(alpha_cruise_3*n^2.*W0S/(0.5*rho_3*V_inf_3max^2*pi*AR*e)));

%%%%%%%%%%%%%%%%%%%%       Absolute ceiling (spec point)        %%%%%%%%%%%%%%%%%%%%%


%40,000ft amsl
sigma_min = 0.2465;
beta_min = 0.225/sigma_min;
Tw0_ceil = (2/(beta_min*sigma_min))*sqrt((Cd0)/(e*pi*AR));


%%%%%%%%%%%%%%%%%%%%       Plotting       %%%%%%%%%%%%%%%%%%%%%
yidx=[1.1,1.1,1.1]; % Y axis block colour
y1=[0.705 0.705];
y2=[1.1 1.1];
x=[-100 15000];
xidx=[5324,7500,15000];


% existing aircraft (An-158, BAE-146, E-190, A318)
wingloadexisting = [4909,5353,5491,5380];
thrustweightexisting = [0.313,0.301,0.35,0.3];

figure(1)
hold on
plot(W0S,TW0_cruise_1,'-','LineWidth',1.5)
plot(W0S,TW0_cruise_2,'-','LineWidth',1.5)
plot(W0S,TW0_cruise_3,'-','LineWidth',1.5)
plot(W0S,TW0_cruise_missed,'-','LineWidth',1.5,'LineWidth',1.5)
area(W0S,TW0_cruise_1max,'FaceColor',[0 0 0.1],'EdgeColor',[0 0 0.1],'LineWidth',1) %0 0.9 0.5
alpha(0.1)
area(W0S,TW0_TO_11,'FaceColor',[0 0.9 1],'EdgeColor',[0 0.9 1],'LineWidth',1)
alpha(0.1)
plot(W0S,TW0_TO_12,'-.','LineWidth',1.5)
plot(W0S,TW0_TO_13,'-.','LineWidth',1.5)
area(xidx,yidx,'FaceColor',[0.8 0.9 0.2],'EdgeColor',[0.8 0.8 0.3],'LineWidth',1)
%plot([W0S_emergency_1,W0S_emergency_1],[0,1],'--')
alpha(0.1)
plot([W0S_11,W0S_11],[0,1],':','LineWidth',2)
plot([W0S_12,W0S_12],[0,1],':','LineWidth',2)
plot([W0S_13,W0S_13],[0,1],':','LineWidth',2)
yline(first_segment_climb,'--','color',[0.4 0.7 0.6],'LineWidth',1.5)
yline(second_segment_climb,'--','color',[0.4 0.3 0],'LineWidth',1.5)
yline(third_segment_climb,'--','color',[1 0.6 0.1],'LineWidth',1.5)
yline(approach_goaround,'--','color',[0.5 0.2 0.4],'LineWidth',1.5)
yline(landing_goaround,'--','color',[0 0.3 0.6],'LineWidth',1.5)
yline(Tw0_ceil,'k--','LineWidth',1.5)%'color',[0.9 0 0],'LineWidth',1.5)
%patch([x fliplr(x)], [y1 fliplr(y2)],[0.8 0 0])
alpha(0.8)
plot(W0S, TW0_cruise_turn,'color',[0.8 0.3 0.6],'LineWidth',1.5)
plot(W0S, TW0_loiter_turn,'color',[0.5 0.8 0.6],'LineWidth',1.5)
p=plot(4850,0.3195,'kd',MarkerFaceColor='k');
%label(p,'Design Point', 'FontSize',10)
alpha(0.1)
%plot(W0S,TW0_cruise_1max,'color','LineWidth',1)
plot([W0S_emergency_1,W0S_emergency_1],[0,1],'color',[0.8 0.8 0.3],'LineWidth',1)


%plot([W0S_emergency_1,W0S_emergency_1],[0,1],'Color',[0.8 0.9 0.2])
%yline(0.705,'color',[0.7 0.7 0.7],'LineWidth',1.5)


p1=plot(4909,0.313,'ko');
label(p1,'An-158', 'FontSize',10)
p2=plot(5353, 0.291,'ko');
label(p2,'BAE-146', 'FontSize',10)
p3=plot(5491, 0.345,'ko');
label(p3,'E-190', 'FontSize',10)
% p4=plot(5380, 0.3,'ko');
% label(p4,'A318')

ylabel('Thrust to Weight (T_o/W_o)', 'FontSize',11)
ylim([0 1])
%title('Aircraft Constraint Diagram', 'FontSize',10)
xlabel('Wing Loading (W_o/S) [N/m^2]', 'FontSize',11)
xlim([0 W0S_11 + 1000])
box on 
grid on

legend('Cruise 1','Cruise 2', 'Cruise 3', 'Missed','Cruise 1 (M = 0.75)', 'Takeoff 1','Takeoff 2', 'Takeoff 3','Emergency landing','Landing 1','Landing 2','Landing 3','1st Segment Climb (OEI)','2nd Segment Climb (OEI)','3rd Segment Climb (OEI)','Go-around Approach (OEI)','Go-around Landing (AEO)', 'Absolute Ceiling', 'Cruise Turn 12.5^o', 'Loiter Turn 25^o','Design Point','Location' , 'eastoutside', 'FontSize',10)
hold off


figure(2)
hold on
plot(W0S,TW0_TO_1,'-.','LineWidth',1.5)
plot(W0S,TW0_TO_11,'-.','LineWidth',1.5)
plot(W0S,TW0_TO_21,'-.','LineWidth',1.5)
plot([W0S_emergency,W0S_emergency],[0,1],'--','LineWidth',1.5)
plot([W0S_emergency_1,W0S_emergency_1],[0,1],'--','LineWidth',1.5)
plot([W0S_emergency_2,W0S_emergency_2],[0,1],'--','LineWidth',1.5)
p=plot(4850,0.3195,'kd',MarkerFaceColor='k');
label(p,'Design Point', 'FontSize',10)
grid on
box on
% yline(TW0_climb,'k')
% yline(Tw0_ceil,'color',[0.9290 0.6940 0.1250])
% yline(0.705,'b')
% %scatter(wingloadexisting,thrustweightexisting)
% ylabel('T_o/W_o', 'FontSize',14)
ylim([0 0.7])


xlabel('W_o/S (N/m^2)', 'FontSize',11)
ylabel('T_o/W_o', 'FontSize',11)
xlim([0 10000])


legend('Takeoff 1 (C_{Lmax}= 2)', 'Takeoff 1 (C_{Lmax}= 2.6)', 'Takeoff 1 (C_{Lmax}= 3)','Emergency landing (C_L_m_a_x= 2)','Emergency landing (C_L_m_a_x= 2.6)','Emergency landing (C_L_m_a_x= 3)','Location' , 'eastoutside', 'FontSize',10)
% 'Cruise 1','Cruise 2', 'Cruise 3', 'Missed','Cruise 1 max', 'Cruise 2
% max', 'Cruise 3 max' 'Climb', 'Absolute Ceiling', 'Max T_o/W_o'
hold off


