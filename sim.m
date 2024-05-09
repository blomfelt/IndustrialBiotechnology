clc
clear
close all

global  q1_max q2_max q3_max q4_max q5_max q6_max q7_max ... 
        K_1Sec K_2G K_2ATP K_3G K_3ATP K_3IATP K_4Pyr ... 
        K_4O2 K_4ISec K_5Pyr K_6E K_6O2 K_6ISec K_7E K_7ATP K_7ISec ...
        gamma21 gamma71 g21 g71 ... 
        kLa P R T He Q ... 
        VL VG yCO2_in yO2_in rho_c 

t       = [0, 85];    % cultivation period (h)
t_Plot  = t;          % To use in the plots later, since ode function changes it

rho_c   = 500;                      % cell density (g DW/(L cell)) 
Vtot    = 100;                      % total volume (L)
VL      = 75;	                    % liquid volume (L)	
VG      = Vtot-VL;                  % gas volume (L)
Q       = 60*VL;                    % gas flow rate (L/h) 

kLa     = 600;                      % oxygen mass transfer coefficient (/h) vart ska vi trycka in kla? 
He      = 0.790;                    % henry's constant for oxygen, (atm*L/mmol) 
P       = 1;                        % pressure (atm)
T       = 293;                      % temperature (K)
R       = 8.206e-5;                 % ideal gas constant (atm*L/mmol*K)
yO2_in  = 0.2095;                   % initial mole fraction of O2
yCO2_in = 0.0005;                   % initial mole fraction of CO2

%pO2 = P*yO2_in; 
%pO2air = 
%DOT = pO2/pO2air 


%Starting conditions
S_ec_i  = 1000;                     % initial extracellular glucose concentration (mM) 
Gi      = 0.1;                      % initial intracellular glucose concentration (mM) 
ATPi    = 1;                        % initial intracellular ATP concentration (mM)
Pyri    = 0.5;                      % Initial intracellular pyruvate concentration (mM)
Xi      = 0.1;                      % Initial biomass concentration (g/L)
cO2_Li  = 100*P*yO2_in/He;          % Initial oxygen concentration in the liquid (mM)
yO2i    = yO2_in;                   % Initial mole fraction of oxygen (1)
yCO2i   = yCO2_in;                  % Initial mole fraction of carbon dioxide (1)
Ei      = 0;                        % Initial ethanol concentration (mM)

%% Kinetic parameters and rate equations
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

starting_conditions = [S_ec_i, Gi, ATPi, Xi, Pyri, cO2_Li, yO2i, yCO2i, Ei]; % initial mole fraction of O2
    
[t, y] = ode15s(@(t,y)sim_fun(t, y), t, starting_conditions);

S_ec   = y(:,1);
G      = y(:,2);
ATP    = y(:,3);
X      = y(:,4);
Pyr    = y(:,5);
CO2_L  = y(:,6);
O2     = y(:,7);
CO2    = y(:,8);
E      = y(:,9);

subplot(3, 2, 1)    
plot(t,X,t,G,'LineWidth',0.75)
title("Concentration of X and G")
legend("biomass, X", "intracellular glucose, G")
ylabel('concentration [g/L]')
xlabel('time [h]')
xlim(t_Plot)

subplot(3, 2, 2)    
plot(t,S_ec,'LineWidth',0.75) 
title("Concentration of S_e_c")
legend("extracellular glucose, S_e_c")
ylabel('concentration [mM]') 
xlabel('time [h]')
xlim(t_Plot)

subplot(3, 2, 3)    
plot(t,ATP,t,Pyr,'LineWidth',0.75)  
title("Concentration of ATP and Pyruvate")
legend('intracellular ATP','intracellular pyruvate')
ylabel('concentration [mM]') 
xlabel('time [h]')
xlim(t_Plot)

subplot(3, 2, 4)    
plot(t,CO2_L,'LineWidth',0.75) 
title("Concentration of CO_2_,_L")
legend('carbon dioxide, CO_2_,_L')
ylabel('concentration [mM]') 
xlabel('time [h]')
xlim(t_Plot)

subplot(3, 2, 5)    
plot(t,E,'LineWidth',0.75)
title("Concentration of E")
legend("ethanol, E")
ylabel('concentration [mM]')
xlabel('time [h]')
xlim(t_Plot)

subplot(3, 2, 6)    
plot(t,O2,'LineWidth',0.75) 
title("Concentration of O_2")
legend("oxygen, O_2")
ylabel('concentration [mM]') 
xlabel('time [h]')
xlim(t_Plot)
