clc
clear
close all

global Fin_l Fin_g V_gas_initial V_initial ...
       yo2_in yco2_in He P T R...
       S0 mu_max kLa Ks ...
       Yxs Yos Ycs Ysx

    t= [0 36]; %h, may also use: linspace(0, 36, 1000);
    
    V_initial=0.5; % L
    V_gas_initial= 9.5; %L
    
    Fin_l=0.2; % L/h
    Fin_g=60;   % L/h
    
    S0= 25; % g/L
    Si=5;   % g/L
    yo2_in=0.2095;
    yco2_in=0.0004;
    
    Xi=0.5;  % g/L
    ci_o2=265*10^-6; % M
    
    mu_max=0.5; % /h
    Ks=0.5; % g/h
    Yxs=0.5; % g/g
    Ysx=1/Yxs;
    Yos=0.01; % mol/g
    Ycs=0.01; % mol/g
    
    kLa=500; % /h
    He=790; % atm*L/mol
    P=1; % atm 
    %P=101325;  % Pa
    T=293;  % K
    R=0.08206; % atm*L/mol*K
    
    
    
    starting_conditions= [Xi, Si, ci_o2, yo2_in, yco2_in];

    [t, y]=ode15s(@(t,y)B4_fun(t, y), t, starting_conditions);
    
    
    X= y(:,1);
    S= y(:,2);
    c_o2= y(:,3);
    y_o2g= y(:,4);
    y_co2g= y(:,5);

subplot(2, 2, 1)
fontScaling = 1.2;
plot(t,X, ...
    t, S)
title("Concentration of X and S")
legend("X", "S")
ylabel("c [g/L]")
%fontsize(gcf, scale=fontScaling)

subplot(2, 2, 2)
plot(t, c_o2*1e6)
ylabel("c(O2) [µ mol]")

subplot(2, 2, 3)
plot(t, y_o2g, t, y_co2g)
title("Molar fractions of O2 and CO2")
legend("y(O2)", "y(CO2)")


V_plot = V_initial+Fin_l*t;
V_gas_plot=V_gas_initial-Fin_l*t;
V_total_plot = V_plot+V_gas_plot;
subplot(2, 2, 4)
plot(t, V_plot, t, V_gas_plot)
ylabel("Volume [L]")
legend("Liquid Volume", "Gas Volume")
fontsize(gcf, scale=fontScaling)
