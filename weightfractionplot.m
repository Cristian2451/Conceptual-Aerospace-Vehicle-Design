clear
clc

A = 0.97;   %constants of regression model
C = -0.06;

Wcrew = 415;    %weights
Wpay = 9085;
Wf_W0 = 0.243;


syms W0     %W0 as symbolic variable to solve eqn

We_W0_reg = A*W0^C; %regression model
We_W0_eqn = - (Wcrew + Wpay)/W0 + 1 - Wf_W0;    %weight eqn

GROSS_TAKEOFF_W = double(vpasolve(We_W0_eqn == We_W0_reg, W0)); %solving intercept

We_W0 = - (Wcrew + Wpay)/GROSS_TAKEOFF_W + 1 - Wf_W0;   %subbing gross weight back into eqn

%plots plots plots
fplot(We_W0_reg,'LineWidth',1.5)
hold on
fplot(We_W0_eqn,'LineWidth',1.5)
hold on
scatter(GROSS_TAKEOFF_W, We_W0, "diamond", "k","filled")
xlim([30000 50000])
xline(GROSS_TAKEOFF_W, '-', ['W_0 = ', num2str(GROSS_TAKEOFF_W)])
yline(We_W0, '-', ['Final W_e/W_0 = ', num2str(We_W0)])
xlabel('Takeoff Weight W_0 (kg)','FontSize',12)
ylabel('Empty Weight Fraction W_e/W_0','FontSize',12)
legend('W_e/W_0 = AW0^c','W_e/W_0 = 1 - W_f/W_0 - W_c/W_0 - W_p/W_0','Location','southeast')
grid on

figure
We = We_W0*GROSS_TAKEOFF_W;
weights = [Wcrew, Wpay, Wf_W0*GROSS_TAKEOFF_W, We];
Wf = Wf_W0 * GROSS_TAKEOFF_W;
pie(weights,{['CREW', newline, num2str(round(Wcrew*100/GROSS_TAKEOFF_W,1)),'%'], ['PAYLOAD', newline, num2str(round(Wpay*100/GROSS_TAKEOFF_W,1)),'%'], ['EMPTY', newline, num2str(round(We*100/GROSS_TAKEOFF_W,1)),'%'], ['FUEL',  newline, num2str(round(Wf*100/GROSS_TAKEOFF_W,2)),'%']})
colormap(gray)
set(gca, 'FontName', 'cmr12')


