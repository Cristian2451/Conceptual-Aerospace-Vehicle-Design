clear
clc

AR = 10;
WSdes = 4850;
W = 39290*9.81;
Sref = W/WSdes;
b = sqrt(AR*Sref);  % wingspan
taper = 0.4;
sweep_c4 = 12;

CLmax_airfoil = 1.55;

CLmax_clean = 0.9*CLmax_airfoil*cosd(sweep_c4);

%% Flaps (Fowler)

dcLmax_base = 1.36; %for 25% chord flap
k1 = 1; %for 25% chord flap
k2 = 1; % 40deg 
k3 = 1;

dcLmf = k1*k2*k3*dcLmax_base;

N = 100;
i = 70; %end of flaps (60% of wing span)
j = 10; %start of flaps (5% of wing span)
theta = pi/(2*N):pi/(2*N):pi/2;
Croot = 2*Sref/((taper+1)*b);
c_start_flap = Croot*(1-(1-taper)*cos(theta(N-j)));
c_end_flap = Croot*(1-(1-taper)*cos(theta(N-i)));

S_flapped = b*((i-j)/N)*(c_start_flap+c_end_flap)/2;
K_lambda = 0.91;

dCLmax_flaps = dcLmf*K_lambda*S_flapped/Sref;


%% Slats

cL_delta_max = 1.2; %for 10% chord ratio
delta_f = pi/12; % 15 deg deflection
Eta_max = 1.56; %for 0.12 LER/(t/c)
Eta_delta = 1; %for delta_f = 15 deg
cdash_c = 1.1; 

dcLms = cL_delta_max*Eta_max*Eta_delta*delta_f*cdash_c;

f = 99; %end of slats (90% of wing span)
g = 10; %start of slats (90% of wing span)
c_start_slat = Croot*(1-(1-taper)*cos(theta(N-g)));
c_end_slat = Croot*(1-(1-taper)*cos(theta(N-f)));

S_slats = b*((f-g)/N)*(c_start_slat+c_end_slat)/2;
sweep_H_L = 13.76; %sweep at hinge line of slats

dCLmax_slats = dcLms*cosd(sweep_H_L)*S_slats/Sref;


%% Overall C_Lmax

CLmax = CLmax_clean + 0.95*dCLmax_flaps + dCLmax_slats;







