clc
clear
close all

global Fin_l Fin_g V_g_i V_l_i ...
       yo2_in yco2_in ci_o2 He P T R...
       S0 mu_max kLa Ks Yxs Yos Ycs Ysx ...
       S_ec G ATP X Pyr CO2 cO2_L E...
       S_eci Gi ATP_i Xi Pyr_i yco2_in yo2_in Ei...
       q1_max q2_max q3_max q4_max q5_max q6_max q7_max...
       K_1Sec g21 g71 gamma71 gamma21 K_2G K_2ATP K_3IATP K_3G K_3ATP...
       K_4Pyr K_4O2 K_4ISec K_5Pyr K_6E   K_6O2  K_6ISec K_7E   K_7ATP K_7ISec

% Kinetic parameters and rate equations

%R1 Glucose uptake
q1_max = 14; %mmol/gDW/h
K_1Sec = 1; %mM

%R2 Growth on glucose
%Stochiometry
g21=10; %mmol ATP/mmol glucose
gamma21=0.15; %gX / mmol glucose, True biomass yield on glucose in biosynthesis
%Kinetics
q2_max = 2.7; %mmol/gDW/h
K_2G = 0.05; %mM
K_2ATP = 0.20; %mM

%R3 Glycolysis
q3_max = 60; %mmol/gDW/h 
K_3G = 0.8; %mM intracell.
K_3ATP = 0.5;%mM intracell.
K_3IATP = 1; %mM intracell.

%R4 Respiration of pyruvate
q4_max = 10; %mmol/gDW/h
K_4Pyr   = 0.2; %mM intracell.
K_4O2   = 0.02; %mM
K_4ISec = 1; %mM 

%5 Fermentation
q5_max = 40; %mmol/gDW/h 
K_5Pyr = 5;  %mM intracell.

%6 Respiration of ethanol and glyoxylate shunt
q6_max = 6; %mmol/gDW/h
K_6E     = 3; %mM
K_6O2   = 0.02; %mM
K_6ISec = 0.5; %mM

%R7 Cell growth on ethanol
%Stoichiometry
g71=12; %mmol ATP/mmol E
gamma71=0.025; %gX/mmol E, True biomass yield on ethanol in biosynthesis
%Kinetics
q7_max = 2; %mmol/gDW/h
K_7E   = 0.5; %mM
K_7ATP = 0.5; %mM
K_7ISec = 0.5; %mM


%S_ec, G,  ATP,       X, Pyr, CO2, O2,  E
%A=[-1, +1,    0,       0,   0,   0,  0,  0;
%    0, -1, -g21, gamma21,   0,   0,  0,  0;
%    0, -1,    2,       0,   2,   0,  0,  0;
%    0,  0,    6,       0,  -1,   3, -3,  0;
%    0,  0,    0,       0,  -1,   1,  0,  1;
%    0,  0,    8,       0,   0,   2, -4, -1;
%    0,  0, -g71, gamma71,   0,   0,  0, -1];

A = [-1,       0,  0,  0,  0,  0,       0;   % S_ec       
      1,      -1, -1,  0,  0,  0,       0;   % G    
      0,    -g21,  2,  6,  0,  8,    -g71;   % ATP    
      0, gamma21,  0,  0,  0,  0, gamma71;   % X    
      0,       0,  2, -1, -1,  0,       0;   % Pyr    
      0,       0,  0,  3,  1,  2,       0;   % CO2    
      0,       0,  0, -3,  0, -4,       0;   % O2
      0,       0,  0,  0,  1, -1,      -1];  % E


S_eci = 1000;         % mmol/L
Gi = 0.1    ;         % mmol/L
ATP_i = 1   ;         % mmol/L
Pyr_i = 0.5 ;         % mmol/L
Xi = 0.1    ;         % g/L
V_l_i = 75  ;         % L
V_g_i = 25  ;         % L
kLa = 600   ;         % /h
He = 0.790  ;         % atm*L/mmmol for oxygen
yo2_in = 1000*(0.2095*100*60)/24.07;      % mmol/h
yco2_in = 1000*(0.0005*100*60)/24.07;      %  mmol/h
R = 8.206e-5;         % atm*L/(mmol*K)
t = linspace(0, 100);  % cultivation period (h)
Ei = 0;

%Fin_l=
%Fin_g=



starting_conditions = [S_eci, Gi, ATP_i, Xi, Pyr_i, yco2_in, yo2_in, Ei];
    
[t,y] = ode15s(@(t,y) sim_fun(t,y), t, starting_conditions);
    
 S_ec  = y(:,1);
 G     = y(:,2);
 ATP   = y(:,3);
 X     = y(:,4);
 Pyr   = y(:,5);
 CO2   = y(:,6);
 cO2_L = y(:,7);
 E     = y(:,8);

subplot(3, 2, 1)    
plot(t,S_ec)
title("S_ec")
legend("S_ec")


subplot(3, 2, 2)    
plot(t,X,t,G) 
title("X and G")
legend("X","G")

subplot(3, 2, 3)    
plot(t,ATP,t,Pyr)  
title("ATP and Pyr")
legend("ATP","Pyr")
ylim([0 10]);


subplot(3, 2, 4)    
plot(t,cO2_L) 
title("cO2_L")
legend("cO2_L")

subplot(3, 2, 5)    
plot(t,CO2) 
title("CO2")
legend("CO2")

subplot(3, 2, 6)    
plot(t,E) 
title("E")
legend("E")