function dydt = sim_fun(t, y)
    global q1_max q2_max q3_max q4_max q5_max q6_max q7_max ... 
           K_1Sec K_2G K_2ATP K_3G K_3ATP K_3IATP K_4Pyr ... 
           K_4O2 K_4ISec K_5Pyr K_6E K_6O2 K_6ISec K_7E K_7ATP K_7ISec ...
           gamma21 gamma71 g21 g71 Vtot... 
           kLa P R T He Q Fin_l Fin_g V_g_i V_l_i... 
           yCO2_in yO2_in rho_c S_ec_i S_in u ...
    
    % Extract state variables from y
    S_ec    = y(1);
    G       = y(2);
    ATP     = y(3);
    X       = y(4);
    Pyr     = y(5);
    cO2_L   = y(6);  
    yO2     = y(7);
    yCO2    = y(8);
    VL      = y(9);
    VG      = y(10);
    E       = y(11);

   
    % Update VL and VG
    %VL = V_l_i + Fin_l * t;  % liquid volume (L)
    %VG = Vtot - VL;          % gas volume (L)
    Fin_g = 60 * VL;             % gas flow rate (L/h)
    kLa = 250 + (75000 - VL) / 300;  % oxygen mass transfer coefficient (/h)

    if cO2_L<0
        cO2_L=0;
    end

    % Kinetic rate expressions
    q1 = q1_max * S_ec / (K_1Sec + S_ec);
    q2 = q2_max * G * ATP / ((K_2G + G) * (K_2ATP + ATP));
    q3 = q3_max * G * ATP / ((K_3G + G) * (K_3ATP + ATP) * (1 + ATP / K_3IATP));
    q4 = q4_max * Pyr * cO2_L / ((K_4Pyr + Pyr) * (K_4O2 + cO2_L) * (1 + S_ec / K_4ISec)); 
    q5 = q5_max * Pyr / (K_5Pyr + Pyr);
    q6 = q6_max * E * cO2_L / ((K_6E + E) * (K_6O2 + cO2_L) * (1 + S_ec / K_6ISec));
    q7 = q7_max * E * ATP / ((K_7E + E) * (K_7ATP + ATP) * (1 + S_ec / K_7ISec));

    % Metabolic rates
    my = (gamma21 * q2 + gamma71 * q7);
    q_G = (q1 - q2 - q3);
    q_ATP = (-g21 * q2 + 2 * q3 + 6 * q4 + 8 * q6 - g71 * q7);
    q_Pyr = (2 * q3 - q4 - q5);
    q_O2 = (-3 * q4 - 4 * q6);
    q_CO2 = (3 * q4 + q5 + 2 * q6);
    q_E = (q5 - q6 - q7);

    % Debugging
    if 0 %t<0.1
        disp('Calculating dS_ecdt');
        disp(['Size of V_l_i: ', num2str(size(V_l_i))]);
        disp(['Size of Fin_l: ', num2str(size(Fin_l))]);
        disp(['Size of VL: ', num2str(size(VL))]);
        disp(['Size of S_in: ', num2str(size(S_in))]);
        disp(['Size of S_ec: ', num2str(size(S_ec))]);
        disp(['Size of q1: ', num2str(size(q1))]);
        disp(['Size of X: ', num2str(size(X))]);
    end
    

     % Differential equations
    dS_ecdt = (Fin_l / VL) * (S_in - S_ec) - q1 * X;
    dGdt = q_G * rho_c - my * G; 
    dATPdt = q_ATP * rho_c - my * ATP;  
    dXdt = my * X;  
    dPyrdt = q_Pyr * rho_c - my * Pyr; 
    dcO2_Ldt = q_O2 * X + kLa * (yO2 * P / He - cO2_L);        % change in oxygen concentration in the liquid phase
    dyO2dt = Fin_g / VG * (yO2_in - yO2 * (1 - yO2_in - yCO2_in) / (1 - yO2 - yCO2)) - kLa * (yO2 * P / He - cO2_L) * VL * R * T / (VG * P);   %  change in mole fraction of oxygen  
    dyCO2dt = Fin_g / VG * (yCO2_in - yCO2 * (1 - yO2_in - yCO2_in) / (1 - yO2 - yCO2)) + q_CO2 * X * VL * R * T / (VG * P);             %  change in mole fraction of carbon dioxide 
    dVLdt = Fin_l;%V_l_i + Fin_l * t;
    dVGdt = -Fin_l;% Vtot - VL;
    dEdt = q_E * X; 

    dydt = [dS_ecdt; dGdt; dATPdt; dXdt; dPyrdt; dcO2_Ldt; dyO2dt; dyCO2dt; dVLdt; dVGdt; dEdt];
end
