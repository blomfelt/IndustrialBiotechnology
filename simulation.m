clc
clear
close all

%% Kinetic parameters and rate equations

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
K_4Pyr = 0.2; %mM intracell.
K_4O2 = 0.02; %mM
K_4ISec = 1; %mM 

%5 Fermentation
q5_max = 40; %mmol/gDW/h 
K_5Pyr = 5;  %mM intracell.

%6 Respiration of ethanol and glyoxylate shunt
q6_max = 6; %mmol/gDW/h
K_6E = 3; %mM
K_6O2 = 0.02; %mM
K_6ISec = 0.5; %mM

%R7 Cell growth on ethanol
%Stoichiometry
g71=12; %mmol ATP/mmol E
gamma71=0.025; %gX/mmol E, True biomass yield on ethanol in biosynthesis
%Kinetics
q7_max = 2; %mmol/gDW/h
K_7E = 0.5; %mM
K_7ATP = 0.5; %mM
K_7ISec = 0.5; %mM

%% funderingar


%X, S_ec, E, O2, CO2, G, ATP, Pyr


R1 = -S_ec*q1 + G*q1;
R2 = -G*q2 - g21*ATP*q2 + gamma21*X*q2
R3 = -G*q3 + 2*Pyr*q3 + 2*ATP*q3
R4 = -3*O2*q4 - Pyr*q4 + 6*ATP*q4 + 3*CO2*q4
R5 = -Pyr*q5 + CO2*q5 + E*q5
R6 = -E*q6 - 4*O2*q6 + 2*CO2*q6 + 8*ATP*q6
R7 = -E*q7 - g71*ATP*q7 + gamma71*X*q7


S_ec = G
G = -g21*ATP + gamma21*X
%% A
%S_ec, G,  ATP,       X, Pyr, CO2, O2,  E
 [-1, +1,    0,       0,   0,   0,  0,  0;
   0, -1, -g21, gamma21,   0,   0,  0,  0;
   0, -1,    2,       0,   2,   0,  0,  0;
   0,  0,    6,       0,  -1,   3, -3,  0;
   0,  0,    0,       0,  -1,   1,  0,  1;
   0,  0,    8,       0,   0,   2, -4, -1;
   0,  0, -g71, gamma71,   0,   0,  0, -1]

%dS_ecdt = -q1
%dGdt = q1-q2-q3
%dATPdt = -g21*q2+2*q3+6*q4+8*q6-g71*q7
%dXdt = gamma21*q2+gamma71*q7
%dPyrdt = 2*q3-q4-q5
%dCO2dt = 3*q4+q5+2*q6
%dO2dt = -3*q4-4*q6
%dEdt = q5-q6-q7

%todo: change these
starting_conditions = [Xi, Si, ci_o2, yo2_in, yco2_in];
    
[t, y]=ode15s(@(t,y)B4_fun(t, y), t, starting_conditions);
      

