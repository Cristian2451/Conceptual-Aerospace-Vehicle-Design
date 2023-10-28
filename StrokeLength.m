clear
clc

g = 9.81;
Vv = 10 * 0.3048;          %m/s vertical descent velocity
eff_shock = 0.9;           %Shock efficiency (~0.65 - 0.9 according to handout)
eff_tyre = 0.47;           %Tyre efficiency according to hand out
Ng = 2.85;                 %Gear load factor (between 2.7 and 3 according to hand out)
MLW = 36884;               %(kg)
W_land = MLW * g;          %Maximum landing weight (N)
P = 1.034e6;               %Pressure in pascals
MTOW = 39290;              %kg

NosePerc = 0.1726; % percentage weight on nose gear statically
H = (0.18 + 2.786/2) + 1.4; %Height of CoG in z above the static floor (m)
B = 15.6471; %Wheel base (m)

%P = 1.034e6;

w = 29.21;                %Tyre width cm
d = 76.2;                 %Tyre diameter cm
W_w = (1 - NosePerc) * MTOW * 9.81 * 1.07 / 4; %Force on each wheel

%Rolling Radius Calculation
Ap = W_w / P;
R_roll = -( ( ( 100^2 * Ap) / (2.3 * sqrt(w * d))) - d / 2) / 100;  %Rolling Radius (m)


R_tyre = d /(2 * 100);           %Radius of the tyre (m)
S_T = R_tyre - R_roll;     %Stroke length of tyre 

%Stroke Length Calculation
S = (Vv^2) / (2*g*eff_shock*Ng) - (S_T * eff_tyre / eff_shock) + (1 * 0.024);   %Required stroke length (m) according to p32 of u/c handout

%Oleo Load Calculation using conversion of energy
L_oleo = (W_land * Vv^2) / (2*g*eff_shock*S);            %From p32 of u/c handout (N)

%Required Oleo Diameter
D_oleo = (0.04 * sqrt(0.5 * L_oleo * 0.2247)) * 0.0254;                %Requires units in lbs but gives diameter in metres

Length_oleo = 2.5 * S;

Length_uc = Length_oleo + R_tyre;

D_piston = D_oleo / 1.3;               % piston diameter in cm


%For Nose Gear in Braking

%W_dyn = (10 * (H * 3.2808) * (MLW * 2.20462) / (B * 3.28084)); %Dynamic load on nose brake (lbs)
W_dyn = 81148 * 0.2248089431;
D_oleo_nose = 0.04 * sqrt(W_dyn) * 0.0254;
D_piston_nose = D_oleo_nose / 1.3;



