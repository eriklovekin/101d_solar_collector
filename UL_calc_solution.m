function UL = UL_calc_solution(Tc1,Tc2,Tpm)
    data = extract_data();
    [n_c,kl,n_ri,L_back,k_back,L_cover,W,len_tube,diam_tube,len_collector,C_b,L_plate,U,eps_req,cp_w,rho_w,mu_w,nu_w,Pr_w,k_w,k_c] = given_vals();

    Ta = data(:,7) + 273; % ambient temp
    V = data(:,21);% windspeed

    Ts = Ta;
    Tp = Tpm;
    Tc2int= Tc2;
    Tc1int= Tc1;
    sigma=5.67e-8;%S-B constant;
    g=9.8;% gravitational const;
    epsilonc=4.*n_ri./(n_ri+1).^2;% emissivity of cover glass; same for covers 1 and 2;
    epsilonp=0.088;% emissiivty of plate. Figure 4.8.3, HW2p1
    d = L_cover;
    Pr = 0.7;% TODO: update with problem-specific value
    nu = 1.96e-5;% kinematic viscosity of air
    k=0.0293;% thermal conductivity of air
    angle = 25.9./180.*pi;%beta in radians
    for m = 1:length(V)
        hw(m,1) = max([5,8.6.*V(m).^0.6./len_collector.^0.4]);% DB 3.15.10
    end
    disp('UL_calc')

    % % From HW3p1, returns UL = Ut = 2.3349 as in hw solution
    % Ts=10+273;% sky T
    % Ta=Ts;% air T = sky T in this problem
    % Tp=100+273; % plate T
    % Tc2int=(Ts+Tp)./2;% initial guess on cover-1 temperature.
    % Tc1int=(Tc2int+Tp)./2;% initial guess on cover-1 temperature.
    % sigma=5.67e-8;%S-B constant;
    % g=9.8;% gravitational const
    % epsilonc=0.8;% emissivity of cover glass; same for covers 1 and 2;
    % epsilonp=0.1;% emissiivty of plate.
    % d=0.025;% gap between c1 AND c2, and between p and c1;
    % Pr=0.7;%Prandtl number;
    % nu=1.96e-5;% kinematic viscosity of air
    % k=0.0293;% thermal conductivity of air
    % angle=45./180.*pi;% tile angle 45 deg, in radiant.
    % hw=10; % enternal convection HTC due to wind
    % k_back = 0;% no back heat loss in hw
    % L_back = 1;%arbitrary since no heat loss

    i=0;
    j=0;
    conv = 1;
while (1)
    Tc2=Tc2int;
    while (1)
        disp(i)
        Tc1=Tc1int;
        hrc2a=sigma.*epsilonc.*(Tc2+Ts).*(Tc2.^2+Ts.^2).*(Tc2-Ts)./(Tc2-Ta);
        hcc2a=hw;
        h1=hcc2a+hrc2a;
        R1=1./h1;

        hrc1c2=sigma.*(Tc1.^4-Tc2.^4)./(Tc1-Tc2)./(1./epsilonc+1./epsilonc-1);
        Rac1c2=g.*(Tc1-Tc2).*1./((Tc1+Tc2)./2).*d.^3.*Pr./nu.^2;
        Nuc1c2=1+max(0,1.44.*(1-1708.*(sin(1.8.*angle)).^1.6./Rac1c2./cos(angle)).*(1-1708./Rac1c2./cos(angle)))+max(0,((Rac1c2.*cos(angle)./5830).^(1./3)-1));
        hcc1c2=Nuc1c2.*k./d;
        h2=hcc1c2+hrc1c2;
        R2=1./h2;

        hrpc1=sigma.*(Tp.^4-Tc1.^4)./(Tp-Tc1)./(1./epsilonp+1./epsilonc-1);
        Rapc1=g.*(Tp-Tc1).*1./((Tp+Tc1)./2).*d.^3.*Pr./nu.^2;
        Nupc1=1+max(0,1.44.*(1-1708.*(sin(1.8.*angle)).^1.6./Rapc1./cos(angle)).*(1-1708./Rapc1./cos(angle)))+max(0,((Rapc1.*cos(angle)./5830).^(1./3)-1));
        hcpc1=Nupc1.*k./d;
        h3=hcpc1+hrpc1;
        R3=1./h3;

        Ut=(R1+R2+R3).^(-1);

        qloss=Ut.*(Tp-Ta);
        Tc1=Tp-qloss./h3;
        Tc2=Tc1-qloss./h2;
        disp([Tc1, Tc2])
        if all(abs(Tc1-Tc1int)<conv)
            Tc1int=Tc1;
            break
        else
            Tc1int=Tc1;% break the while loop if Tc1 is close to Tc1int
            i=i+1;
        end
    end
    if all(abs(Tc2-Tc2int)<conv)
        Tc2int=Tc2;
        break% break the while loop if Tc2 is close to Tc2int
    else
        Tc2int=Tc2;
        j=j+1;
    end

    Ub = k_back./L_back;
    UL = real(Ut + Ub);
end