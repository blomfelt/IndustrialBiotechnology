clc
clear
close all

global  Fin_l Fin_g V_g_i V_l_i ...
        yo2_in yco2_in ci_o2 He P T R ... 
        S0 mu_max kLa Ks  ...
        Yxs Yos Ycs Ysx 

    t = linspace(0, 36, 1000); % cultivation period (h)

    V_l_i   = 0.5;             % initial liquid volume (L)
    V_g_i   = 9.5;             % initial gas volume (L)
    
    Fin_l   = 0.02;             % liquid flow rate (L/h)
    Fin_g   = 60;              % gas flow rate (L/h) 

    VL      = V_l_i+Fin_l*t;   % liquid volume, at time t (L) 
    VG      = V_g_i-Fin_l*t;   % gas volume, at time t (L)

    S0      = 100;              % substrate concentration in inlet (g/L)
    Si      = 5;               % initial substrate concentration (g/L)
    yo2_in  = 0.2095;          % mole fraction of O2 in inlet
    yco2_in = 0.0004;          % mole fraction of CO2 in inlet 
    
    Xi      = 0.5;             % initial biomass concentration (g/L)
    ci_o2   = 265e-6;          % initial O2 concentration, liquid (M) 
    
    mu_max  = 0.5;             % maximum specific growth rate (/h)
    Ks      = 0.5;             % substrate saturation constant (g/L)
    Yxs     = 0.5;             % biomass yield (g biomass/g substrate)
    Ysx     = 1/Yxs;           % substrate yield (g substrate/ g biomass)
    Yos     = 0.01;            % oxygen yield (mol/g) 
    Ycs     = 0.01;            % carbon dioxide yield (mol/g)
    
    kLa     = 500;             % oxygen mass transfer coefficient (/h)
    He      = 790;             % henry's constant for oxygen, (atm*L/mol) 
    P       = 1;               % pressure (atm)
    T       = 293;             % temperature (K)
    R       = 0.08206;         % ideal gas constant (atm*L/mol*K)
    
    starting_conditions = [Xi, Si, ci_o2, yo2_in, yco2_in];
    
    [t, y]=ode15s(@(t,y)B4_fun(t, y), t, starting_conditions);
      
    X= y(:,1);
    S= y(:,2);
    c_o2= y(:,3);
    y_o2g= y(:,4);
    y_co2g= y(:,5);

subplot(2, 2, 1)    
plot(t,X,t,S,'LineWidth',0.75)
title("Concentration of X and S")
legend("biomass concentration, X", "substrate concentration, S")
ylabel('concentration [g/L]')
xlabel('time [h]')
xlim([0 36])

subplot(2, 2, 2)    
plot(t,c_o2*10^6,'LineWidth',0.75) % ändra skala, är det mikro nu?
title("Oxygen concentration in the liquid")
legend("c_O_2_,_l")
ylabel('concentration [μM]') 
xlabel('time [h]')
xlim([0 36])

subplot(2, 2, 3)    
plot(t,y_o2g,t,y_co2g,'LineWidth',0.75)  
title("Mole fractions of oxygen and carbon dioxide in outlet")
legend('y_O_2_,_g','y_C_O_2_,_g')
ylabel('mole fraction') 
xlabel('time [h]')
xlim([0 36])

subplot(2, 2, 4)    
plot(t,VL,t,VG,'LineWidth',0.75) 
title("Volume")
legend('liquid','gas')
ylabel('volume [L]') 
xlabel('time [h]')
xlim([0 36])

