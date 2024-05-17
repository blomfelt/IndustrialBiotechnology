clc
clear
close all

global q1_max q2_max q3_max q4_max q5_max q6_max q7_max ... 
       K_1Sec K_2G K_2ATP K_3G K_3ATP K_3IATP K_4Pyr ... 
       K_4O2 K_4ISec K_5Pyr K_6E K_6O2 K_6ISec K_7E K_7ATP K_7ISec ...
       gamma21 gamma71 g21 g71 Vtot... 
       P R T He Fin_l Fin_g ... 
       yCO2_in yO2_in rho_c S_ec_i S_in u V_l_i V_g_i

tmax    = 100;         % total hours
%tspan   = [0 tmax];    % cultivation period (h)
dt      = 0.5;
tspan   = 0:dt:tmax;

rho_c   = 500;                      % cell density (g DW/(L cell)) 

V_l_i   = 20000;             % initial liquid volume (L)
V_g_i   = 80000;             % initial gas volume (L)

%Fin_l   = 4000;             % liquid flow rate (L/h)
Fin_l    = 0; % Start at zero

%Fin_g   = 60;              % gas flow rate (L/h) 
Vtot    = 100000;           % total volume (L)
VL_max  = 75000;            % maximum liquid volume (L)

He      = 0.790;                    % henry's constant for oxygen, (atm*L/mmol) 
P       = 1;                        % pressure (atm)
T       = 293;                      % temperature (K)
R       = 8.206e-5;                 % ideal gas constant (atm*L/mmol*K)
yO2_in  = 0.2095;                   % initial mole fraction of O2
yCO2_in = 0.0005;                   % initial mole fraction of CO2

% Starting conditions
S_in    = 1000;                      % glucose concentration (mmol/L) inflow 
S_ec_i  = 500;                       % initial extracellular glucose concentration (mmol/L) 
Gi      = 0.1;                      % initial intracellular glucose concentration (mM) 
ATPi    = 1;                        % initial intracellular ATP concentration (mM)
Pyri    = 0.5;                      % Initial intracellular pyruvate concentration (mM)
Xi      = 0.375;                    % Initial biomass concentration (g/L)
cO2_Li  = 100 * P * yO2_in / He; 
yO2i    = yO2_in; 
yCO2i   = yCO2_in; 
Ei      = 0;                        % Initial ethanol concentration ()

% Kinetic parameters and rate equations
%R1 Glucose uptake
q1_max  = 14;                       % (mmol/gDW/h)
K_1Sec  = 1;                        % (mM)
%R2 Growth on glucose
%Stochiometry
g21     = 10;                       % (mmol ATP/mmol glucose)
gamma21 = 0.15;                     % (gX / mmol glucose) True biomass yield on glucose in biosynthesis
%Kinetics
q2_max  = 2.7;                      % (mmol/gDW/h)
K_2G    = 0.05;                     % (mM)
K_2ATP  = 0.20;                     % (mM)
%R3 Glycolysis
q3_max  = 60;                       % (mmol/gDW/h)
K_3G    = 0.8;                      % (mM) intracell
K_3ATP  = 0.5;                      % (mM) intracell
K_3IATP = 1;                        % (mM) intracell
%R4 Respiration of pyruvate
q4_max  = 10;                       % (mmol/gDW/h)
K_4Pyr  = 0.2;                      % (mM) intracell
K_4O2   = 0.02;                     % (mM)
K_4ISec = 1;                        % (mM) 
%5 Fermentation
q5_max  = 40;                       % (mmol/gDW/h)
K_5Pyr  = 5;                        % (mM) intracell
%6 Respiration of ethanol and glyoxylate shunt
q6_max  = 6;                        % (mmol/gDW/h)
K_6E    = 3;                        % (mM)
K_6O2   = 0.02;                     % (mM)
K_6ISec = 0.5;                      % (mM)
%R7 Cell growth on ethanol
%Stoichiometry
g71     = 12;                       % (mmol ATP/mmol E)
gamma71 = 0.025;                    % (gX/mmol E) True biomass yield on ethanol in biosynthesis
%Kinetics
q7_max  = 2;                        % (mmol/gDW/h)
K_7E    = 0.5;                      % (mM)
K_7ATP  = 0.5;                      % (mM)
K_7ISec = 0.5;                      % (mM) 

starting_conditions = [S_ec_i, Gi, ATPi, Xi, Pyri, cO2_Li, yO2i, yCO2i, V_l_i, V_g_i, Ei]; 

% PID control parameters
%Kp = 0.025;
%Ki = 0.016;
%Kd = 0.009;
Kp = 0.15;

E_setpoint = 5;     % setpoint ethanol concentration

% Avoiding errors
y = zeros(length(tspan), 11);
E_error = zeros(length(tspan), 1);
y(1, :) = starting_conditions;
E_error(1) = E_setpoint;
PID_active = 0;
u=[];

for k = 2:length(tspan)
    [t, y_out] = ode15s(@sim_fun, [tspan(k-1), tspan(k)], y(k-1, :));
    y(k, :) = y_out(end, :);

    % Turn on PID:
    if y(k,1) < S_ec_i && y(k,11) > E_setpoint/2 % Turn on the PID control once glucose and ethanol has been produced
%    if y(k,1) < S_ec_i % Turn on the PID control glucose has been consumed
%    if 1
    PID_active = 1; 
    end

    % PID
    if PID_active
        E_error(k) = E_setpoint - y(k, 11); % error between setpoint and current ethanol concentration

        proportional = E_error(k); % calculate proportional
        %integral = sum(E_error);   % calculate integral
        %derivative = (E_error(k-1) - E_error(k)) / dt;  % calculate derivative
        u(k) = Kp * proportional;% + Ki * integral + Kd * derivative;  % Calculate the control signal
        Fin_l = Fin_l - u(k);             
    else
        Fin_l = 0;
    end
    
    if Fin_l<0
        disp("Avoiding negative flow in")
        Fin_l = 0;
    end
    if y(k, 9) > VL_max 
        disp("Max volume reached")
        y(k, 9) = VL_max;
        y(k, 10) = Vtot-VL_max;
        Fin_l=0;
        Fin_g=0;
    end

% Debugging
  if k<20
      disp('Fin_l:'); disp(Fin_l);
      %disp('u:'); disp(u);
  end
end

% Extract results
S_ec   = y(:, 1);
G      = y(:, 2);
ATP    = y(:, 3);
X      = y(:, 4);
Pyr    = y(:, 5);
cO2_L  = y(:, 6);  
yO2    = y(:, 7);
yCO2   = y(:, 8);
VL     = y(:, 9);
VG     = y(:, 10);
E      = y(:, 11);




t_plot = linspace(0, tmax, length(y));
% Plotting results
subplot(3, 2, 1)    
plot(t_plot, X, t_plot, G, 'LineWidth', 0.75)
title("Concentration of X and G")
legend("biomass, X", "intracellular glucose, G")
ylabel('concentration [g/L]')
xlabel('time [h]')
xlim([0 tmax])

subplot(3, 2, 2)    
plot(t_plot, S_ec, 'LineWidth', 0.75) 
title("Concentration of S_e_c")
legend("extracellular glucose, S_e_c")
ylabel('concentration [mmol/L]') 
xlabel('time [h]')
xlim([0 tmax])

subplot(3, 2, 3)    
plot(t_plot, ATP, t_plot, Pyr, 'LineWidth', 0.75)  
title("Concentration of ATP and Pyruvate")
legend('intracellular ATP','intracellular pyruvate')
ylabel('concentration [mM]') 
xlabel('time [h]')
xlim([0 tmax])

subplot(3, 2, 4)    
plot(t_plot, cO2_L * 10^3, 'LineWidth', 0.75) 
title("O_2_,_L")
legend('O_2_,_L')
ylabel('enhet]') 
xlabel('time [h]')
xlim([0 tmax])

subplot(3, 2, 5)    
plot(t_plot, E, 'LineWidth', 0.75)
title("Concentration of E")
legend("ethanol, E")
ylabel('concentration enhet')
xlabel('time [h]')
xlim([0 tmax])


subplot(3, 2, 6)    
plot(t_plot, yO2, t_plot, yCO2) 
title(" O_2 and CO_2")
legend("O_2","CO_2")
ylabel('concentration enhet') 
xlabel('time [h]')
ylim([0 1])
xlim([0 tmax])

figure()
subplot(2, 2, 1)
plot(1:length(u), u)
title("u")
xlabel("Timepoint (not hour!)")

subplot(2, 2, 2)
plot(1:length(E_error), E_error)
title("E error")
xlabel("Timepoint (not hour!)")

subplot(2, 2, 3)
plot(1:length(VL), VL)
title("VL")
xlabel("Timepoint (not hour!)")

