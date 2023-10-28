clear
clc

%% Takeoff Distance

T_W = 0.33;
W_S = 4850;
AR = 10;
Wo = 39290*9.81;
miu = 0.03; 
CD0_takeoff = 0.098;
CL_Takeoff = 1.96; 
e = 0.764; 
rho = 1.225; 

V_stall = 56.85;
V_LOF = 1.1*V_stall;

K_T = T_W - miu;
K_A = (rho/(2*W_S))*(miu*CL_Takeoff - CD0_takeoff - (CL_Takeoff^2/(pi*AR*e)));

S_G = (1/(2*9.81*K_A))*log((K_T + K_A*V_LOF^2)/K_T);

S_R = 3*V_LOF; %for a 3 second rotation time (Raymer)


V_TR = 1.15*V_stall;
n = 1.2;   %load factor (from Raymer)
L_D = 0.5*sqrt(pi*AR*e/CD0_takeoff);  

R = V_TR^2/((n-1)*9.81);
gamma_cl = asin(T_W - 1/L_D);

hobs = 10.7;
htr = R*(1-cos(gamma_cl));

if hobs >= htr
    S_TR = R*sin(gamma_cl);
    S_CL = (hobs - htr)/tan(gamma_cl);
else
    S_TR = sqrt(R^2 - (R - hobs)^2);
    S_CL = 0;
end

S_TO = 1.15*(S_G + S_R + S_CL + S_TR); %Total takeoff distance


%% Balanced Field Length

CL_climb = 1.96; 
sigma = 1; 
CLmax_clean = 1.365; 
BPR = 4.9;
Ne = 2;

U = 0.01*CLmax_clean + 0.02;
D_OEI = 0; 
gamma_min = 0.024;
gamma_cl_OEI = 1.8;%asin(((Ne - 1)/Ne)*T_W - D_OEI/Wo);
G = gamma_cl_OEI - gamma_min;
T_mean = 0.75*64540*((5 + BPR)/(4 + BPR));
T_W_mean = T_mean/Wo;

BFL = (0.863/(1+2.3*G))*(W_S/(rho*9.81*CL_climb + hobs))*(1/(T_W_mean - U) + 2.7) + 655/sqrt(sigma);


%% Landing Distance

miu_landing = 0.5;
CL_landing = 0; 
CD0_landing = 0.2499;
%L_D_approach = 0.5*sqrt(pi*AR*e/CD0_landing); 

V_a = 1.3*V_stall;
V_F = 1.23*V_stall;
V_TD = 1.15*V_stall;

T_approach = 0.15*64540*2;

gamma_a = 3*pi/180;%asin(T_approach/Wo - 1/L_D_approach);

R = V_F^2/((n-1)*9.81);

hf = R*(1 - cos(gamma_a));

S_A = (hobs - hf)/tan(gamma_a);

S_F = R*sin(gamma_a);

S_FR = 3*V_TD;

T_landing = (0.15 -0.4)*64540*2; 

K_T = T_landing/Wo - miu_landing;
K_A = (rho/(2*W_S))*(miu_landing*CL_landing - CD0_landing - (CL_landing^2/(pi*AR*e)));

S_B = (1/(2*9.81*K_A))*log((K_T)/(K_T + K_A*V_TD^2));

S_L = 1.666*(S_A + S_F + S_FR + S_B);

% Without thrust reversers

T_landing_NTR = (0.15)*64540*2; 

K_T_NTR = T_landing_NTR/Wo - miu_landing;
K_A_NTR = (rho/(2*W_S))*(miu_landing*CL_landing - CD0_landing - (CL_landing^2/(pi*AR*e)));

S_B_NTR = (1/(2*9.81*K_A_NTR))*log((K_T_NTR)/(K_T_NTR + K_A_NTR*V_TD^2));

S_L_NTR = 1.666*(S_A + S_F + S_FR + S_B_NTR);



%% Mission Performance

CD0 = 0.0223;
MTOW = 39290;
MF = 9587;
MZFW = MTOW - MF;
V_cruise1 = 143*sqrt(1/0.4481);
V_cruise2 = 136*sqrt(1/0.3475);
V_cruise3 = 127.6*sqrt(1/0.395);
V_cruiseDiv = 124*sqrt(1/0.3881);
V_loiter = 96*sqrt(1/0.8617);

R_cruise1 = 360; % km
R_cruise2 = 710;
R_cruise3 = 420;
R_cruiseDiv = 370;
R_loiter = 45*60*V_loiter/1000;

W_END_LAND = 29797;
L_D = 0.5*sqrt(pi*AR*e/CD0);
sfc = 0.47/3600;
R_extra = -(200/sfc)*L_D*log(MZFW/W_END_LAND)/1000;

payload = 9975;

R_MAXFUEL_LESSPAYLOAD = -(200/sfc)*L_D*log((MZFW-800)/MZFW)/1000;

R_NO_PAYLOAD = -(200/sfc)*L_D*log((MZFW-payload)/MZFW)/1000;

R_MZFW = R_cruise1 + R_cruise2 + R_cruise3 + R_cruiseDiv + R_loiter + R_extra;
R_MAX = R_MZFW + R_NO_PAYLOAD;

FUELSAVED_W_NO_PAYLOAD = 2505.03;

R_FUELSAVED_W_NO_PAYLOAD = -(200/sfc)*L_D*log((MZFW-payload)/(MZFW+FUELSAVED_W_NO_PAYLOAD-payload))/1000;


R_800KGEXTRAFUEL = -(200/sfc)*L_D*log((MZFW-800)/MZFW)/1000;


figure(1)
xline(R_cruise1,'b-.',LineWidth=1)
hold on
xline(R_cruise1+R_cruise2,'r-.',LineWidth=1)
xline(R_cruise1+R_cruise2+R_cruise3,'g-.',LineWidth=1)
xline(R_cruise1+R_cruise2+R_cruise3+R_cruiseDiv,'c-.',LineWidth=1)
xline(R_cruise1+R_cruise2+R_cruise3+R_cruiseDiv+R_loiter,'k-.',LineWidth=1)
plot([R_MZFW R_MZFW+R_FUELSAVED_W_NO_PAYLOAD],[payload 0],'--',LineWidth=1,Color=[0 0.4470 0.7410])
plot([R_MZFW R_MAX],[payload 0],':',LineWidth=1.5,Color=[0 0.4470 0.7410])
plot([0 R_MZFW],[payload payload],'-',LineWidth=1.5,Color=[0 0.4470 0.7410])
plot([R_MZFW R_MZFW+R_800KGEXTRAFUEL+100],[payload payload-800],'-',LineWidth=1.5,Color=[0 0.4470 0.7410])
plot([R_MZFW+R_800KGEXTRAFUEL+100 R_MZFW+R_FUELSAVED_W_NO_PAYLOAD+R_800KGEXTRAFUEL+100],[payload-800 0],'-',LineWidth=1.5,Color=[0 0.4470 0.7410])
legend('End Cruise 1','End Cruise 2','End Cruise 3','End Cruise Diversion','End Loiter','Reduced Payload','Reduced Payload $\rightarrow$ Increased Fuel','Interpreter', 'latex', 'FontSize',12)
xlabel('Range [km]','Interpreter', 'latex', 'FontSize',14)
ylabel('Payload [kg]','Interpreter', 'latex', 'FontSize',14)
ylim([0 10500])
xlim([0 13000])
grid on
box on


%% Flight Envelope

h = linspace(0,13000,100);
[T, a, P, rho] = atmosisa(h);

CD0 = 0.025; 
K = 1/(pi*AR*e);
Mcrit = 0.69;

% At low subsonic speeds

K1 = 1;
K2 = 0;
K3 = -0.6;
K4 = -0.04;
BPR = 4.9;


V_v = linspace(0,4000,9)*0.00508;


syms V

% Stall line
CL_max = 2.45;
V_s = sqrt(2*W_S./(rho*CL_max));
M_s = V_s./a;
figure(2)
plot(M_s,h*3.28084,'r--',LineWidth=1)
hold on
xline(0.741,'k-.',LineWidth=1)
yline(40000,'b:',LineWidth=2)

for j = 1:9
    for i = 1:100
        if h(i) < 11000
            s = 0.7;
        else
            s = 0.7;
        end

        T_To = (K1 + K2*BPR + (K3 + K4*BPR)*(V/a(i)))*(rho(i)/1.225)^s;

        eqn = V*(T_W*T_To - 0.5*rho(i)*V^2*CD0/W_S - n^2*K*W_S/(0.5*rho(i)*V^2)) == V_v(j);

        Velocity = double(solve(eqn,V));
        
        M(:,i,j) = Velocity(:)/a(i);
        
        for k = 1:length(Velocity)
            if M(k,i,j) > 0.385
                M(k,i,j) = NaN;
            end 
        end
        if abs(M(2,i,j) - M(3,i,j)) < 0.1
            M(2,i,j) = NaN;
            M(3,i,j) = NaN;
        end
    end
    plot(M(:,:,j),h*3.28084,LineWidth=1,Color=[0 0.4470 0.7410])
end

clear M;
clear V;
clear Velocity;
clear T_To;
clear eqn;

% At high subsonic speeds

K1 = 0.88;
K2 = -0.016;
K3 = -0.3;
K4 = 0;
s = 0.7;
BPR = 4.9;

syms V


for j = 1:9
    for i = 1:100
        if h(i) < 11000
            s = 0.7;
        else
            s = 0.7;
        end
        T_To = (K1 + K2*BPR + (K3 + K4*BPR)*(V/a(i)))*(rho(i)/1.225)^s;

        eqn = V*(T_W*T_To - 0.5*rho(i)*V^2*CD0/W_S - n^2*K*W_S/(0.5*rho(i)*V^2)) == V_v(j);

        Velocity = double(solve(eqn));
        M(:,i,j) = Velocity(:)/a(i);

        if M(3,i,j) > Mcrit
            %CD0 = 0.021 + 20*(M(2,i,j) - Mcrit)^4;
            T_To = (K1 + K2*BPR + (K3 + K4*BPR)*(V/a(i)))*(rho(i)/1.225)^s;
            eqn = V*(T_W*T_To - 0.5*rho(i)*V^2*CD0/W_S - n^2*K*W_S/(0.5*rho(i)*V^2)) == V_v(j);
            Velocity = double(solve(eqn,V));
            M(2,i,j) = Velocity(3)/a(i);
        end

        for k = 1:length(Velocity)
            if M(k,i,j) < 0.415
                M(k,i,j) = NaN;
            end 
        end
        out = true;
        if i > 1
            out = isreal(M(2,i-1,j));
        end
        
        if out == false
            M(2,i,j) = NaN;
            M(3,i,j) = NaN;
            
        end
    end

    p1 = plot(M(2,:,j),h*3.28084,LineWidth=1,Color=[0 0.4470 0.7410]);
    p2 = plot(M(3,:,j),h*3.28084,LineWidth=1,Color=[0 0.4470 0.7410]);
    if j == 1
        label(p1,'0')
        %label(p2,'0')
    elseif j == 2
        label(p1,'500')
        %label(p2,'500')
    elseif j == 3
        label(p1,'1000')
        %label(p2,'1000')
    elseif j == 4
        label(p1,'1500')
        %label(p2,'1500')
    elseif j == 5
        label(p1,'2000')
        %label(p2,'2000')
    elseif j == 6
        label(p1,'2500')
        %label(p2,'2500')
    elseif j == 7
        label(p1,'3000')
        %label(p2,'3000')
    elseif j == 8
        label(p1,'3500')
        %label(p2,'3500')
    elseif j == 9
        label(p1,'4000')
        %label(p2,'4000')
    end

end
xlim([0 1])
ylim([0 42500])
ylabel('Altitude [ft]','Interpreter', 'latex', 'FontSize',14)
xlabel('Mach','Interpreter', 'latex', 'FontSize',14)
legend('Stall line','Drag divergence Mach number','Absolute ceiling requirement','Vertical velocity lines [ft/min]','Interpreter', 'latex', 'FontSize',12)
grid on
box on
hold off

