clear
clc


%LOWER CASE L IN CL IS FOR SECTION
%UPPER CASE L IN CL IS FOR WING

W_S = 4850; %wing loading
W_W0 = [0.96 0.869 0.766 0.722];    %weight fractions at start of each cruise
Vi = [143 136 127.6 124];     %IAS for each cruise
rho = 1.225;    %sea level density
Mcruise = 0.69;     %cruise Mach

CL_req = (2*W_S.*W_W0./(rho.*Vi.^2)).*sqrt(1-Mcruise^2); %required CL for each leg
%CL_req_avg = (CL_req(1) + CL_req(2) + CL_req(3) + CL_req(4))/4; %average of all cruises


Mcrit12 = @(Cl) -0.1299*Cl + 0.723; %Mcrit models for NACA64012
%Mcrit10 = @(Cl) -0.1299*Cl + 0.741;    %NACA64010



x_torenbeek = 0.5;  %correction exponent for cosine rule

phi = 20;   %initial guess
error = 100;
its = 0;    %initialise iterations counter

while abs(error)>0.00001
    Cl_req = CL_req(1)/(0.9*cosd(phi)^2);    % calculating section Cl with estimated 10% losses due to 3d effects
    %tc_eff = 0.12 * cosd(phi);  %effective t/c based on cosine rule
    %Mcrit_aerofoil = ((0.12 - tc_eff)/(0.12 - 0.1))*(Mcrit10(Cl_req) - Mcrit12(Cl_req)) + Mcrit12(Cl_req);  %linear interp of Mcrit for t/c eff
    Mcrit_aerofoil = Mcrit12(Cl_req);
    phi_ = acosd((Mcrit_aerofoil/Mcruise)^(1/x_torenbeek)); %updating sweep angle based on corrected cosine rule
    error = phi - phi_; %update error
    phi = phi_; %assign new phi to old phi
    its = its+1;    %update iterations counter
end

phi
