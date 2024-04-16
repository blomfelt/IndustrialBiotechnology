function dydt = B4_fun(t, y, VG, VL, Fin_l, Fin_g, S0, V_gas_initial, V_initial,yco2_in, mu_max, kLa, Ks, yo2_in, He, Ysx, Yxs, Yos, Ycs, P, T, R)
  
%global t, y, Fin_l, Fin_g, V_gas_initial, V_initial, yco2_in, S0, mu_max, kLa, Ks, yo2_in, He, Yxs, Yos, Ycs, Ysx, P, T, R

    X = y(1);
    S = y(2);
    c_o2 = y(3);
    yo2 = y(4);
    yco2 = y(5);


    mu=mu_max*S/(Ks+S);
    q_o2=kLa*(yo2_in-yo2)./(VL*X);
    q_co2=-kLa*(yo2_in-yo2)./(VL*X);

    dV_Ldt = Fin_l;          % change in liquid volume
    dV_Gdt = -Fin_l;         % change in gas volume
    dXdt = -Fin_l*X+mu*X;    % change in X
    dsdt =  Fin_l./VL*(S0-S)-Ysx*mu*X;      % substrate consumption rate
    dc_o2dt = q_o2*X+kLa*((yo2*P)./He-c_o2);     % O2 in liquid phase
    dyo2dt =  Fin_g./VG.*(yo2_in-yo2*(1-yo2_in-yco2_in)/(1-yo2-yco2))-kLa*(yo2*P/He-c_o2)*VL*R*T./(VG*P);   %  O2 in gas 
    dyco2dt = Fin_g./VG.*(yco2_in-yco2*(1-yo2_in-yco2_in)/(1-yo2-yco2))+q_co2*X.*VL*R*T./(VG*P);    %  CO2 i gas

   %dV_Ldt; dV_Gdt;
    dydt=[dXdt; dsdt; dc_o2dt; dyo2dt; dyco2dt];
    