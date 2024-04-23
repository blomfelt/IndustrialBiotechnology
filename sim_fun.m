function dydt = sim_fun(t, y)
%% old
global Fin_l Fin_g V_g_i V_l_i ...
       yo2_in yco2_in ci_o2 He P T R...
       S0 mu_max kLa Ks ...
       Yxs Yos Ycs Ysx

%todo: change these to the correct ones
    X       = y(1);
    S       = y(2);
    c_o2    = y(3);
    yo2     = y(4);
    yco2    = y(5);

    VL      = V_l_i+Fin_l*t;
    VG      = V_g_i-Fin_l*t;

    if VL>10
        VL = 10; 
        VG = 0; 
        Fin_l=0;
        Fin_g=0;
    end

    %% new
%Kinetic rate expressions:
q1=q1_max*S_ec/((K_1Sec+S_ec));
q2=q2_max*G*ATP/((K_2G+G)*(K_2ATP+ATP));
q3=q3_max*G*ATP/((K_3G+G)*(K_3ATP+ATP)*(1+ATP/K_3IATP));
q4=q4_max*Pyr*cO2_L/((K_4Pyr+Pyr)*(K_4O2+cO2_L)*(1+S_ec/K_4ISec)); 
q5=q5_max*Pyr/(K_5Pyr+Pyr);
q6=q6_max*E*cO2_L/((K_6E+E)*(K_6O2+cO2_L)*(1+S_ec/K_6ISec));
q7=q7_max*E*ATP/((K_7E+E)*(K_7ATP+ATP)*(1+S_ec/K_7ISec));


dS_ecdt = -q1
dGdt = q1-q2-q3
dATPdt = -g21*q2+2*q3+6*q4+8*q6-g71*q7
dXdt = gamma21*q2+gamma71*q7
dPyrdt = 2*q3-q4-q5
dCO2dt = 3*q4+q5+2*q6
dO2dt = -3*q4-4*q6
dEdt = q5-q6-q7
