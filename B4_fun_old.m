function dydt = B42_fun(t, y)

global Fin_l Fin_g V_gas_initial V_initial ...
       yo2_in yco2_in He P T R...
       S0 mu_max kLa Ks ...
       Yxs Yos Ycs Ysx

    X = y(1);
    S = y(2);
    c_o2 = y(3);
    yo2 = y(4);
    yco2 = y(5);


    VL=V_initial+Fin_l*t;
    VG=V_gas_initial-Fin_l*t;

if VL>10
    VL = 10;
    VG = 0;
    Fin_l=0;
    Fin_g=0;
end


    mu=mu_max*S/(Ks+S);
    q_o2= -Yos/Yxs*mu; %kLa*(yo2_in-yo2)/(VL*X); %OLD
    q_co2= Ycs/Yxs*mu;

    dV_Ldt = Fin_l;          % change in liquid volume
    dV_Gdt = -Fin_l;         % change in gas volume
    dXdt = -Fin_l/VL*X+mu*X;    % change in X
    dsdt =  Fin_l/VL*(S0-S)-Ysx*mu*X;      % substrate consumption rate
    dc_o2dt = q_o2*X+kLa*((yo2*P)/He-c_o2);     % O2 in liquid phase
    dyo2dt =  Fin_g/VG*(yo2_in-yo2*(1-yo2_in-yco2_in)/(1-yo2-yco2))-kLa*(yo2*P/He-c_o2)*VL*R*T/(VG*P);   %  O2 in gas 
    dyco2dt = Fin_g/VG*(yco2_in-yco2*(1-yo2_in-yco2_in)/(1-yo2-yco2))+q_co2*X*VL*R*T/(VG*P);    %  CO2 i gas
    
    dydt=[dXdt; dsdt; dc_o2dt; dyo2dt; dyco2dt];
    