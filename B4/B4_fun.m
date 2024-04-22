function dydt = B4_fun(t, y)

global Fin_l Fin_g V_g_i V_l_i ...
       yo2_in yco2_in ci_o2 He P T R...
       S0 mu_max kLa Ks ...
       Yxs Yos Ycs Ysx

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

    mu      = mu_max*S/(Ks+S);                   % specific growth rate (/h)
    q_o2    = -Yos/Yxs*mu;                       % specific rate of oxygen consumption (mol/g*h)
    q_co2   = Ycs/Yxs*mu;                        % specific rate of carbon dioxide production (mol/g*h)   

    dV_Ldt  = Fin_l;                             % change in liquid volume
    dV_Gdt  = -Fin_l;                            % change in gas volume 
    dXdt    = -Fin_l/VL*X+mu*X;                  % change in biomass concentration
    dsdt    = Fin_l/VL*(S0-S)-Ysx*mu*X;          % change in substrate concentration 
    dc_o2dt = q_o2*X+kLa*(yo2*P/He-c_o2);        % change in oxygen cocnentration in the liquid phase
    dyo2dt  = Fin_g/VG*(yo2_in-yo2*(1-yo2_in-yco2_in)/(1-yo2-yco2))-kLa*(yo2*P/He-c_o2)*VL*R*T/(VG*P);   %  change in mole fraction of oxygen  
    dyco2dt = Fin_g/VG*(yco2_in-yco2*(1-yo2_in-yco2_in)/(1-yo2-yco2))+q_co2*X*VL*R*T/(VG*P);             %  change in mole fraction of carbon dioxide
    dydt    = [dXdt; dsdt; dc_o2dt; dyo2dt; dyco2dt];
    