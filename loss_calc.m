%% problem 1

L = 0.025; % plate to cover spacing
ep= 0.1; % plate emittance
Ta = 283.15; % ambient temp
hw = 10; % wind ht coefficient
% coltilt = 45; % collector tilt in deg
Tfi = 273.15+55; % fluid inlet temp
ec = 0.8; % glass cover emittance


% initial guesses:
Tp = 373.15; % mean plate temp
Tc1 = 350; % cover 1
Tc2 = 300; % cover 2
Ut = 7; 
Ut_new = 10;

% use this guess to find HT coeffs

while abs(Ut_new-Ut)>0.01
    Ut = Ut_new;
    hrpc1=htcoeff(Tp,Tc1,ep,ec);
    Ra = 9.81*0.00367*(Tp-Tc1)*L^3/15.06e-6/21.7e-6;
    Nu = 1+1.44*(1-1708*(sind(1.8*coltilt)^1.6/Ra/cosd(coltilt)))*(1-1708/Ra/cosd(coltilt))+((Ra*cosd(coltilt)/5830)^(1/3)-1);
    hcpc1 = Nu*0.02587/L;
    
    R3= 1/(hrpc1+hcpc1);
    
    hrc1c2=htcoeff(Tc1,Tc2,ec,ec);
    Ra = 9.81*0.00367*(Tc1-Tc2)*L^3/15.06e-6/21.7e-6;
    Nu = 1+1.44*(1-1708*(sind(1.8*coltilt)^1.6/Ra/cosd(coltilt)))*(1-1708/Ra/cosd(coltilt))+((Ra*cosd(coltilt)/5830)^(1/3)-1);
    hcc1c2 = Nu*0.02587/L;
    
    R2 = 1/(hrc1c2+hcc1c2);
    
    hrc2a=5.67e-8*(Tc2+Ta)*(Tc2^2+Ta^2)*(Tc2-Ta)/(Tc2-Ta);
    
    R1 = 1/(hrc2a+hw);
    
    m=sqrt((Ut_new+0.9)/385/5e-4)^0.5;
    F=tanh(m*(0.15-0.01)/2)/m/(0.15-.01)*2;
    F_prime= 1/(Ut_new+0.9)/.15/(1/(Ut_new+0.9)/(0.01+(0.15-0.01)*F)+1/pi/.01/300);
   

    dimflow = 0.03*4190/2/(Ut_new+0.9)/F_prime;
    F_doubleprime = dimflow*(1-exp(-1/dimflow));
    F_R = F_prime*F_doubleprime;


    ambienttemps = [-11,-8,-2,2,3,6,7,8,9,7];
    tempdiffs = 40-ambienttemps;
    Svalues = [.01,.35,.82,3.29,2.84,3.39,3.21,1.63,0.99,0.04];
    q_u = zeros(1,length(ambienttemps));

    for n=1:length(ambienttemps)
        q_u(n) = F_R*(Svalues(n)*10^6-tempdiffs(n)*(Ut_new+0.9)*3600);
    end

    sumintensity = 19.79e6;
    eff_day = sum(q_u)/sumintensity;
    dailyusefulgain = 10*2*sum(q_u);

    Tp = Tfi + Qu/Ac/F_R/U_L*(1-F_R)
    Tc1 = Tp - Ut*(Tp-Ta)/(hcpc1+hrpc1);
    Tc2 = Tc1 - Ut*(Tp-Ta)/(hcc1c2+hrc1c2);

    % this feels like a wildly inefficient way to do this loop

    hrpc1=htcoeff(Tp,Tc1,ep,ec);
    Ra = 9.81*0.00367*(Tp-Tc1)*L^3/15.06e-6/21.7e-6;
    Nu = 1+1.44*(1-1708*(sind(1.8*coltilt)^1.6/Ra/cosd(coltilt)))*(1-1708/Ra/cosd(coltilt))+((Ra*cosd(coltilt)/5830)^(1/3)-1);
    hcpc1 = Nu*0.02587/L;
    
    R3= 1/(hrpc1+hcpc1);
    
    hrc1c2=htcoeff(Tc1,Tc2,ec,ec);
    Ra = 9.81*0.00367*(Tc1-Tc2)*L^3/15.06e-6/21.7e-6;
    Nu = 1+1.44*(1-1708*(sind(1.8*coltilt)^1.6/Ra/cosd(coltilt)))*(1-1708/Ra/cosd(coltilt))+((Ra*cosd(coltilt)/5830)^(1/3)-1);
    hcc1c2 = Nu*0.02587/L;
    
    R2 = 1/(hrc1c2+hcc1c2);
    
    hrc2a=5.67e-8*(Tc2+Ta)*(Tc2^2+Ta^2)*(Tc2-Ta)/(Tc2-Ta);
    
    R1 = 1/(hrc2a+hw);
    
    Ut_new = 1/(R1+R2+R3);
end

fprintf("1. The top loss coefficient for this solar collector is: "+round(Ut_new,3)+" W/m^2K")