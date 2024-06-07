function UL = UL_calc(Tc1,Tc2,Tpm,Ac)

data = extract_data();

n_covers = 2;
n_tubes = 2; % starting with 2
KL_covers = .0125;
RI_covers = 1.526; % refractive index
back_insulation_thickness = .007;
back_insulation_conductivity = .0245;
spacing = 0.025;
tube_distance = 0.115;
tube_length = 2.46;
tube_diameter = 0.015; % thin walled
collector_length = 2.5;
bond_conductance = 10e6;
plate_thickness = 0.0005;

%% problem 1

L = spacing; % plate to cover spacing
ep= 0.096; % average plate emittance from hw 2
Ta = data(:,7)'; % ambient temp
windspeeds = data(:,end);
Re = 1.232.*windspeeds*collector_length/1.794e-5; % page 339
% from 3.15.11/12:
if Re(:)<1000
    Nu = .4+.54.*Re.^.52;
elseif Re(:)>1000
    Nu = 0.30.*Re.^.6;
else
    fprintf("damn.")
end
hw = KL_covers/collector_length^2.*Nu; % wind ht coefficient
coltilt = 25.9; % collector tilt in deg, found in S calcs
ec = 4*RI_covers/(RI_covers+1)^2; % glass cover emittance, find a credible source for this

Tpm = 373.15; % guessing mean plate temp
Tc1 = 350; % initial guess
Tc2 = 300; %initial guess
Ut = 7; % initial guess
Ut_new = 10;

% use this guess to find HT coeffs

while abs(Ut_new-Ut)>0.01
    Ut = Ut_new;
    hrpc1=htcoeff(Tpm,Tc1,ep,ec);
    Ra = 9.81*0.00367*(Tpm-Tc1)*L^3/15.06e-6/21.7e-6;
    Nu = 1+1.44.*(1-1708.*(sind(1.8*coltilt).^1.6./Ra./cosd(coltilt))).*(1-1708./Ra./cosd(coltilt))+((Ra.*cosd(coltilt)/5830).^(1/3)-1);
    hcpc1 = Nu*0.02587/L;

    R3= 1./(hrpc1+hcpc1);

    hrc1c2=htcoeff(Tc1,Tc2,ec,ec);
    Ra = 9.81*0.00367.*(Tc1-Tc2).*L^3/15.06e-6/21.7e-6;
    Nu = 1+1.44.*(1-1708*(sind(1.8*coltilt).^1.6./Ra/cosd(coltilt))).*(1-1708./Ra./cosd(coltilt))+((Ra.*cosd(coltilt)/5830).^(1/3)-1);
    hcc1c2 = Nu*0.02587/L;

    R2 = 1./(hrc1c2+hcc1c2);

    hrc2a=5.67e-8.*(Tc2+Ta).*(Tc2.^2+Ta.^2).*(Tc2-Ta)./(Tc2-Ta);

    R1 = 1./(hrc2a+hw);

    Tc1 = Tpm - Ut.*(Tpm-Ta)./(hcpc1+hrpc1);
    Tc2 = Tc1 - Ut.*(Tpm-Ta)./(hcc1c2+hrc1c2);

    % this feels like a wildly inefficient way to do this loop

    hrpc1=htcoeff(Tpm,Tc1,ep,ec);
    Ra = 9.81*0.00367*(Tpm-Tc1)*L^3/15.06e-6/21.7e-6;
    Nu = 1+1.44.*(1-1708.*(sind(1.8.*coltilt)^1.6./Ra./cosd(coltilt))).*(1-1708./Ra./cosd(coltilt))+((Ra.*cosd(coltilt)/5830).^(1/3)-1);
    hcpc1 = Nu*0.02587/L;

    R3= 1./(hrpc1+hcpc1);

    hrc1c2=htcoeff(Tc1,Tc2,ec,ec);
    Ra = 9.81*0.00367.*(Tc1-Tc2)*L.^3/15.06e-6/21.7e-6;
    Nu = 1+1.44.*(1-1708.*(sind(1.8.*coltilt).^1.6./Ra./cosd(coltilt))).*(1-1708./Ra/cosd(coltilt))+((Ra.*cosd(coltilt)/5830).^(1/3)-1);
    hcc1c2 = Nu*0.02587/L;

    R2 = 1./(hrc1c2+hcc1c2);

    hrc2a=5.67e-8.*(Tc2+Ta).*(Tc2.^2+Ta.^2).*(Tc2-Ta)./(Tc2-Ta);

    R1 = 1./(hrc2a+hw);

    Ut_new = 1./(R1+R2+R3);
end

    function htc = htcoeff(T1,T2,eps1, eps2)
        htc = (5.67e-8.*(T2.^2+T1.^2).*(T2+T1))./((1-eps1)./eps1+1+(1-eps2)./eps2);
    end

Ut_back = back_insulation_thickness/back_insulation_conductivity % from 6.4.10, k/L

return real(Ut_new+Ut_back)

end