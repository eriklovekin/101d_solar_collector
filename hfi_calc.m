%% Calculate hfi
function hfi = hfi_calc(m_dot_sc,T_sc_out,T_sc_in)
    [n_c,kl,n_ri,L_back,k_back,L_cover,W,len_tube,diam_tube,len_collector,C_b,L_plate,U,eps_req,cp_w,rho_w,mu_w,nu_w,Pr_w,k_w,k_c] = given_vals();
        % Get water values from temp, for now assume constant
        nu = nu_w;
        Pr = Pr_w;
        rho = rho_w;
        cp = cp_w;
        k = k_w;

        D = diam_tube;
        L = len_tube;

        Ac = pi*diam_tube^2/4;%pipe area
        G = m_dot_sc./Ac;%mass velocity BHT 4.39
        Re = G.*D./mu_w;% BHT 4.39
        Re = real(Re);
        
        %check laminar/turbulent %Lec 13-16 notes
        if Re <= 2300 %laminar
            % disp('hfi_calc: laminar')
            % Nu = 3.66 + (0.065*(D/L)*Re*Pr)/(1+0.04*((D/L)*Re*Pr)^(2/3)); %BHT 4.50 averaged for entrance effect
            Nu = 4.364;% DB 4.41
            % f = 64/Re;% BHT 4.39
        else %turbulent
            % disp('hfi_calc: turbulent')
            if Pr < 0.5
                disp("hfi_calc: Nu calc may not be valid because Pr < 0.5")
            end
            % f = (0.79*log(Re)-1.64)^-2;% DB 4.42
            Nu = 0.023.*Re.^0.8.*Pr.*0.4;%BHT 4.44, could also use 4.45 or 3.14.1 from DB
            
        end
        Nu = real(Nu);
        hfi = Nu*k/D;
end
