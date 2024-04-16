clc
clear

%global  Fin_l, Fin_g, V_gas_initial, V_initial, yco2_in, S0, mu_max, kLa, Ks, yo2_in, He, Yxs, Yos, Ycs, Ysx, P, T, R
 t= linspace(0, 36, 1000); % h
    V_initial=0.5; % L
    V_gas_initial= 9.5; %L
    
    Fin_l=0.2; % L/h
    Fin_g=60;   % L/h

    VL=V_initial+Fin_l*t;
    VG=V_gas_initial+Fin_g*t;
    

    
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
    

    
    [t, y]=ode15s(@(t,y)B4_fun(t, y, VG, VL, Fin_l, Fin_g, V_gas_initial, V_initial, yco2_in, S0, mu_max, kLa, Ks, yo2_in, He, Yxs, Yos, Ycs, Ysx, P, T, R), t, starting_conditions);
    
    
    X= y(:,1);
    S= y(:,2);
    c_o2= y(:,3);
    y_o2g= y(:,4);
    y_co2g= y(:,5);
    
    
    plot(t,X)
    
    

