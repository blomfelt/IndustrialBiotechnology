function dydt = sim_fun(t, y)

global q1_max q2_max q3_max q4_max q5_max q6_max q7_max ... 
       K_1Sec K_2G K_2ATP K_3G K_3ATP K_3IATP K_4Pyr ... 
       K_4O2 K_4ISec K_5Pyr K_6E K_6O2 K_6ISec K_7E K_7ATP K_7ISec ...
       gamma21 gamma71 g21 g71 ... 
       kLa P R T He Q ... 
       VL VG yCO2_in yO2_in rho_c 

    S_ec    = y(1);
    G       = y(2);
    ATP     = y(3);
    X       = y(4);
    Pyr     = y(5);
    cO2_L   = y(6);  
    yO2     = y(7);
    yCO2    = y(8);
    E       = y(9);

 %Kinetic rate expressions
    q1      = q1_max*S_ec/((K_1Sec+S_ec));
    q2      = q2_max*G*ATP/((K_2G+G)*(K_2ATP+ATP));
    q3      = q3_max*G*ATP/((K_3G+G)*(K_3ATP+ATP)*(1+ATP/K_3IATP));
    q4      = q4_max*Pyr*cO2_L/((K_4Pyr+Pyr)*(K_4O2+cO2_L)*(1+S_ec/K_4ISec)); 
    q5      = q5_max*Pyr/(K_5Pyr+Pyr);
    q6      = q6_max*E*cO2_L/((K_6E+E)*(K_6O2+cO2_L)*(1+S_ec/K_6ISec));
    q7      = q7_max*E*ATP/((K_7E+E)*(K_7ATP+ATP)*(1+S_ec/K_7ISec));


    my      = (gamma21*q2+gamma71*q7);
    q_G     = (q1-q2-q3);
    q_ATP   = (-g21*q2+2*q3+6*q4+8*q6-g71*q7);
    q_Pyr   = (2*q3-q4-q5);
    q_O2    = (-3*q4-4*q6);
    q_CO2   = (3*q4+q5+2*q6);
    q_E     = (q5-q6-q7);
    

    dS_ecdt     = -q1*X;
    dGdt        = q_G*rho_c - my*G;
    dATPdt      = q_ATP*rho_c - my*ATP;
    dXdt        = my*X;
    dPyrdt      = q_Pyr*rho_c - my*Pyr;
    dcO2_Ldt    = q_O2*X + kLa*(yO2*P/He-cO2_L); % change in oxygen concentration in the liquid phase
    dyO2dt      = Q/VG*(yO2_in-yO2*(1-yO2_in-yCO2_in)/(1-yO2-yCO2))-kLa*(yO2*P/He-cO2_L)*VL*R*T/(VG*P);   %  change in mole fraction of oxygen  
    dyCO2dt     = Q/VG*(yCO2_in-yCO2*(1-yO2_in-yCO2_in)/(1-yO2-yCO2))+q_CO2*X*VL*R*T/(VG*P);             %  change in mole fraction of carbon dioxide
    dEdt        = q_E*X;  

    dydt    = [dS_ecdt; dGdt; dATPdt; dXdt; dPyrdt; dcO2_Ldt; dyO2dt; dyCO2dt; dEdt];
