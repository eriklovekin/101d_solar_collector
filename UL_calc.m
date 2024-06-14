function UL = UL_calc(Tc1,Tc2,Tpm)

data = extract_data();
[n_c,kl,n_ri,L_back,k_back,L_cover,W,len_tube,diam_tube,len_collector,C_b,L_plate,U,eps_req,cp_w,rho_w,mu_w,nu_w,Pr_w,k_w,k_c] = given_vals();

% n_covers = 2;
% n_tubes = 2; % starting with 2
% KL_covers = .0125;
% n_ri = 1.526; % refractive index
% back_insulation_thickness = .007;
% back_insulation_conductivity = .0245;
% spacing = 0.025;
% tube_distance = 0.115;
% tube_length = 2.46;
% tube_diameter = 0.015; % thin walled
% len_collector = 2.5;
% bond_conductance = 10e6;
% plate_thickness = 0.0005;

% rho_air = 1.225;%kg/m^3 Source: https://www.aerodynamics4students.com/properties-of-the-atmosphere/sea-level-conditions.php
% mu_air = 1.789e-5;% kg/m/s
% 
% L = L_cover; % plate to cover spacing
% ep= 0.096; % average plate emittance from hw 2
% Ta = data(:,7) + 273; % ambient temp
% V = data(:,21); %windspeed
% Re = rho_air.*V.*len_collector/mu_air; % DB page 330
% % from 3.15.11/12:
% if Re(:)<1000
%     Nu = .4+.54.*Re.^.52;
% else% Re(:)>1000
%     Nu = 0.30.*Re.^.6;
% end
% % hw = kl/len_collector^2.*Nu; % wind ht coefficient
% for k = 1:length(V)
%     hw(k,1) = max([5,8.6.*V(k).^0.6/len_collector.^0.4]);% DB 3.15.10
% end
% beta = 25.9; % collector tilt in deg, found in S calcs
% ec = 4*n_ri/(n_ri+1)^2; % glass cover emittance, TODO: find a credible source for this

% %Check: hw3 p1 values {valid!}
L = 0.025; % plate to cover spacing
ep= 0.1; % plate emittance
Ta = 283.15; % ambient temp
hw = 10; % wind ht coefficient
Tpm = 373.15; % mean plate temp
beta = 45; % collector tilt in deg
ec = 0.8; % glass cover emittance
Tc1 = 350; % initial guess
Tc2 = 300; %initial guess

Ut = 7; % initial guess
Ut_new = 10;

% use this guess to find HT coeffs

while abs(Ut_new-Ut)>0.01
    % disp('UL_calc')
    Ut = Ut_new;
    hrpc1=htcoeff(Tpm,Tc1,ep,ec);
    Ra = 9.81.*0.00367.*(Tpm-Tc1).*L^3./15.06e-6./21.7e-6;
    Nu = 1+1.44.*(1-1708.*(sind(1.8*beta).^1.6./Ra./cosd(beta))).*(1-1708./Ra./cosd(beta))+((Ra.*cosd(beta)/5830).^(1/3)-1);%%% makes inf
    hcpc1 = Nu*0.02587/L;

    R3= 1./(hrpc1+hcpc1);

    hrc1c2=htcoeff(Tc1,Tc2,ec,ec);
    Ra = 9.81*0.00367.*(Tc1-Tc2).*L^3/15.06e-6/21.7e-6;
    Nu = 1+1.44.*(1-1708*(sind(1.8*beta).^1.6./Ra/cosd(beta))).*(1-1708./Ra./cosd(beta))+((Ra.*cosd(beta)/5830).^(1/3)-1);
    hcc1c2 = Nu*0.02587/L;

    R2 = 1./(hrc1c2+hcc1c2);

    hrc2a=5.67e-8.*(Tc2+Ta).*(Tc2.^2+Ta.^2).*(Tc2-Ta)./(Tc2-Ta);

    R1 = 1./(hrc2a+hw);

    Tc1 = Tpm - Ut.*(Tpm-Ta)./(hcpc1+hrpc1);
    Tc2 = Tc1 - Ut.*(Tpm-Ta)./(hcc1c2+hrc1c2);

    % this feels like a wildly inefficient way to do this loop

    hrpc1=htcoeff(Tpm,Tc1,ep,ec);
    Ra = 9.81*0.00367*(Tpm-Tc1)*L^3/15.06e-6/21.7e-6;
    Nu = 1+1.44.*(1-1708.*(sind(1.8.*beta)^1.6./Ra./cosd(beta))).*(1-1708./Ra./cosd(beta))+((Ra.*cosd(beta)/5830).^(1/3)-1);
    hcpc1 = Nu*0.02587/L;

    R3= 1./(hrpc1+hcpc1);

    hrc1c2=htcoeff(Tc1,Tc2,ec,ec);
    Ra = 9.81*0.00367.*(Tc1-Tc2)*L.^3/15.06e-6/21.7e-6;
    Nu = 1+1.44.*(1-1708.*(sind(1.8.*beta).^1.6./Ra./cosd(beta))).*(1-1708./Ra/cosd(beta))+((Ra.*cosd(beta)/5830).^(1/3)-1);
    hcc1c2 = Nu*0.02587/L;

    R2 = 1./(hrc1c2+hcc1c2);

    hrc2a=5.67e-8.*(Tc2+Ta).*(Tc2.^2+Ta.^2).*(Tc2-Ta)./(Tc2-Ta);

    R1 = 1./(hrc2a+hw);

    Ut_new = 1./(R1+R2+R3);
end

    function htc = htcoeff(T1,T2,eps1, eps2)
        htc = (5.67e-8.*(T2.^2+T1.^2).*(T2+T1))./((1-eps1)./eps1+1+(1-eps2)./eps2);
    end

% Ut_back = L_back/k_back; % from 6.4.10, k/L
Ut_back = k_back/L_back;
UL = real(Ut_new+Ut_back);
UL = Ut_new';% hw1 check

end
