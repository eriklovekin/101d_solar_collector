% Main Loop
clear all
close all

converged = 0.1;
Tci = 20 + 273;% Heat exchanger cold side in temp
Tco = 55 + 273;% Heat exchanger cold side out temp

T_sc_out = (61.2 + 273)*ones([23,1]); % Set Solar Collector OUT!!!

volume_min = 120*0.00378541;%m^3

% Year,Month,Day,Hour,Minute,DHI,Temperature,Clearsky DHI,Clearsky DNI,Clearsky GHI,Cloud Type,Dew Point,DNI,GHI,Relative Humidity,Solar Zenith Angle,Surface Albedo,Pressure,Precipitable Water,Wind Direction,Wind Speed
%   1    2    3   4     5     6       7           8             9           10          11        12     13   14        15                 16              17           18           19               20            21 
[n_c,kl,n_ri,L_back,k_back,L_cover,W,len_tube,D,len_collector,Cb,delta,U,eps_req,cp_w,rho_w,mu_w,nu_w,Pr_w,k_w,k_c] = given_vals();
data = extract_data();
wind_speed = data(:,21);
wind_dir = data(:,20);
Ta = data(:,7) + 273;

n = 1;% number of tubes
increase_area = true;

S = s_calc();

while increase_area
    % Increase area
    n = n+1;
    % disp(['n = ' num2str(n)])

    Ac = len_collector*n*W;

    eps = eps_req*ones(23,1);% Reset effectiveness requirement to eps_req
       
    % Assuming Ch = Cmin, DB 3.17.5 simplifies
    T_sc_in = T_sc_out - eps.*(T_sc_out-Tci);

    %Guess values to calculate m_dot_h
    m_dot_sc = 0.005*ones([23,1]);%kg/s
    m_dot_sc_old = 1;
    % Tpm = (T_sc_out+T_sc_in)/2;% mean plate temp
    Tpm = 336;
    Tpm_old = Tpm + 1;
    % Tc2 = (Ta+Tpm)./2;
    % Tc1 = (Tc2+Tpm)./2;
    Tc1 = 309;
    Tc2 = 298;
    i = 0;

    % disp('   [Rc       Ntu1        Ntu      mdot_he      A       maxA     maxA_loc    eps12  min_eps]') %Legend for mdot_and_A_he display
    %% Based on initial guesses, calculate Tpm and m_dot in solar collector
    while any(abs(m_dot_sc-m_dot_sc_old) > converged) || any(abs(Tpm-Tpm_old) > converged)
        i = i + 1;
        % disp(i);
        
        % UL = UL_calc(Tc1,Tc2,Tpm);
        UL = UL_calc_solution(Tc1,Tc2,Tpm);

        hfi = hfi_calc(m_dot_sc);% convective heat transfer coefficient inside tube (doesn't use Thi,Tho currently)

        m = sqrt(UL./(k_c*delta));% DB 6.5.4a
        F = tanh(m.*(W-D)./2)./(m.*(W-D)./2);% DB 6.5.12 standard fin efficiency
        
        tmp1 = 1./(UL.*(D + (W - D).*F));
        tmp2 = 1./(pi.*D.*hfi);
        F_prime = (1./UL)./(W.*(tmp1 + 1/Cb + tmp2));% DB 6.5.18

        FR = m_dot_sc.*cp_w.*(T_sc_out-T_sc_in)./(Ac.*(S - UL.*(T_sc_in - Ta)));% DB 6.7.1
        for k = 1:length(FR)
            if FR(k) < 0
                FR(k) = 0;
            end
        end
        %update m_dot_h and Tpm
        m_dot_sc_old = m_dot_sc;
        Tpm_old = Tpm;

        tmp3 = (T_sc_out - Ta - S./UL)./(T_sc_in - Ta - S./UL);
        m_dot_sc = (-UL.*Ac.*F_prime)./(log(tmp3)*cp_w);% DB 6.6.4

        m_dot_sc = real(m_dot_sc);% clean up m_dot_sc
        for k = 1:length(m_dot_sc)
            if m_dot_sc(k) < 0
                m_dot_sc(k) = 0;
            end
            if isnan(m_dot_sc(k))
                m_dot_sc(k) = 0;
            end
        end
        
        % Qu = Ac.*(S - UL.*(Tpm - Ta));% DB 6.9.3
        Qu = Ac.*FR.*(S - UL.*(T_sc_in - Ta));% 6.7.6
        Tpm = T_sc_in + (Qu/Ac).*(1-FR)./(FR.*UL);% DB 6.9.4
        
        Tpm = real(Tpm);
        for k = 1:length(Tpm) % clean up Tpm
            if Tpm(k) < 0
                Tpm(k) = 1;
            end
            if isnan(Tpm(k))
                Tpm(k) = 2;
            end
        end

        [m_dot_he, Ahe, eps] = mdot_and_A_he(eps,m_dot_sc,T_sc_out);
        for k = 1:length(eps)
            if isnan(eps(k))
                eps(k) = 1;
            end
        end
        % disp('     UL        hfi       m         F        F_prime    FR       Qu       m_dot_he   m_dot_sc    Ahe       eps       Tpm')
        % disp([UL,hfi,m,F,F_prime,FR,Qu,m_dot_he,m_dot_sc,Ahe*ones(size(eps)),eps,Tpm])
        % disp(['Ahe = ' num2str(Ahe)])
        % Check that initial assumption is met
        if m_dot_he < m_dot_sc
            disp('heat exchanger cold m_dot is less than hot m_dot. Increasing SC area')
            break %increase solar colector area
        end

        %% Using updated eps values, recalc T_sc_in and guess new cover temps
        % Assuming Ch = Cmin, DB 3.17.5 simplifies
        T_sc_in = T_sc_out - eps.*(T_sc_out-Tci);
        Tc2 = (Ta+Tpm)./2;
        Tc1 = (Tc2+Tpm)./2;
    end

    heated_volume = sum(m_dot_he./rho_w)*3600;
    % v_hourly(:,n-8) = m_dot_he./rho_w*3600*1000;
    % v_tot(n-8) = heated_volume*1000;
    hv_gal = heated_volume/0.00378541;

    % disp(['Hot water produced: ', num2str(heated_volume), ' m^3 = ', num2str(hv_gal), ' gallons'])
    if heated_volume >= volume_min
        increase_area = false;
    end
end
disp(['Hot water produced: ', num2str(heated_volume), ' m^3 = ', num2str(hv_gal), ' gallons'])
disp(['Area of collector: ' num2str(Ac)])
disp(['Number of tubes: ' num2str(n)])





