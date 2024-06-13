function [m_dot_he,maxA,eps] = mdot_and_A_he(eps,m_dot_sc)
%% Calculate m_dot_he (cold side of HE), Area of heat exchanger, and update eps values for this new area
    [n_c,kl,n_ri,L_back,k_back,L_cover,W,len_tube,diam_tube,len_collector,C_b,L_plate,U,eps_req,cp_w,rho_w,mu_w,nu_w,Pr_w,k_w,k_c] = given_vals();
    T_sc_out = 61.3 + 273;
    T_he_in = 20 + 273;
    T_he_out = 55 + 273;

    m_dot_he = m_dot_sc.*eps.*(T_sc_out-T_he_in)./(T_he_out-T_he_in); % DB 3.17.5

    Rc = m_dot_sc./m_dot_he; % BHT 8.34
    for k = 1:length(Rc)
            if isnan(Rc(k))
                Rc(k) = 1;
            end
    end
    Ntu1 = 1./(1-Rc).*log((1-eps.*Rc)./(1-eps)); % BHT Table 8.3b, counterflow
    Ntu1 = real(Ntu1);
    for k = 1:length(Ntu1) % cleanup Ntu1
        if isnan(Ntu1(k))
            Ntu1(k) = 0;
        end
        if isinf(Ntu1(k))
            Ntu1(k) = 50;
        end
    end
    A=Ntu1.*m_dot_sc.*cp_w./U; % DB 3.17.8;
    for k = 1:length(A)% cleanup A
        if isnan(A(k))
            A(k) = 0;
        end
    end
    [maxA, max_location]=max(A); % TODO: should always be at noon, but isn't in i = 4 and a few others.
    Ntu = U*A(12)./(m_dot_sc.*cp_w);% DB 3.17.8,  Hardcoded to be at noon
    Ntu = real(Ntu); %cleanup Ntu
    for k = 1:length(Ntu)
        if isnan(Ntu(k))
        Ntu(k) = 0;
        end
        if isinf(Ntu(k))
            Ntu(k) = 50;
        end
    end
    % eps = (1-exp((1-Rc).*Ntu))./(Rc-exp((1-Rc).*Ntu));% DB 3.17.8, Rc = C*, and BHT Table 8.3b        OLD AND MAYBE BAD?
    eps = (1-exp(-Ntu.*(1-Rc)))./(1-Rc.*exp(-Ntu.*(1-Rc)));% DB 3.17.7a, Rc = C*
    eps = real(eps);% cleanup eps
    for k = 1:length(eps)
        if isnan(eps(k))
        eps(k) = 1;
        end
    end
    % disp([Rc(12),Ntu1(12),Ntu(12),m_dot_he(12),A(12),maxA,max_location,eps(12),min(eps)])
end