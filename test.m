n = 1; % number of tubes
increase_area = true;

while increase_area
    n = n + 1;
    disp(['n = ' num2str(n)])
    Ac = len_collector * n * W;

    eps = eps_req * ones(23,1); % Reset effectiveness requirement to eps_req

    % Assuming Ch = Cmin, DB 3.17.5 simplifies
    T_sc_in = T_sc_out - eps .* (T_sc_out - Tci);

    % Guess values to calculate m_dot_h
    m_dot_sc = 0.1 * ones([23,1]); % kg/s
    m_dot_sc_old = 1;
    Tpm = (T_sc_out + T_sc_in) / 2; % mean plate temp
    Tpm_old = Tpm + 1;
    Tc2 = (Ta + Tpm) / 2;
    Tc1 = (Tc2 + Tpm) / 2;
    i = 0;

    % Calculation loop for m_dot_sc and Tpm
    while any(abs(m_dot_sc - m_dot_sc_old) > converged) || any(abs(Tpm - Tpm_old) > converged)
        i = i + 1;
        disp(i);

        UL = UL_calc_solution(Tc1, Tc2, Tpm);
        hfi = hfi_calc(m_dot_sc, T_sc_out, T_sc_in); % Convective heat transfer coefficient inside tube

        m = sqrt(UL / (k_c * delta)); % DB 6.5.4a
        F = tanh(m * (W - D) / 2) / (m * (W - D) / 2); % DB 6.5.12 standard fin efficiency

        tmp1 = 1 / (UL * (D + (W - D) * F));
        tmp2 = 1 / (pi * D * hfi);
        F_prime = (1 / UL) / (W * (tmp1 + 1 / Cb + tmp2)); % DB 6.5.18

        FR = m_dot_sc .* cp_w .* (T_sc_out - T_sc_in) ./ (Ac .* (S - UL .* (T_sc_in - Ta))); % DB 6.7.1

        % Avoid negative FR values
        FR(FR < 0) = 0;

        % Update m_dot_sc and Tpm
        m_dot_sc_old = m_dot_sc;
        Tpm_old = Tpm;

        tmp3 = (T_sc_out - Ta - S ./ UL) ./ (T_sc_in - Ta - S ./ UL);
        m_dot_sc = (-UL * Ac * F_prime) ./ (log(tmp3) * cp_w); % DB 6.6.4

        % Clean up m_dot_sc
        m_dot_sc = real(m_dot_sc);
        m_dot_sc(m_dot_sc < 0) = 0;
        m_dot_sc(isnan(m_dot_sc)) = 0;

        Qu = Ac .* FR .* (S - UL .* (T_sc_in - Ta)); % 6.7.6
        Tpm = T_sc_in + (Qu / Ac) .* (1 - FR) ./ (FR .* UL); % DB 6.9.4

        % Clean up Tpm
        Tpm = real(Tpm);
        Tpm(Tpm < 0) = 1;
        Tpm(isnan(Tpm)) = 2;

        [m_dot_he, Ahe, eps] = mdot_and_A_he(eps, m_dot_sc);
        eps(isnan(eps)) = 1;

        disp('     m_dot_he   eps     Tpm')
        disp([m_dot_he, eps, Tpm])

        % Check if initial assumption is met
        if m_dot_he < m_dot_sc
            disp('heat exchanger cold m_dot is less than hot m_dot. Increasing SC area')
            break % Increase solar collector area
        end

        % Recalculate T_sc_in and cover temps using updated eps values
        T_sc_in = T_sc_out - eps .* (T_sc_out - Tci);
        Tc2 = (Ta + Tpm) / 2;
        Tc1 = (Tc2 + Tpm) / 2;
    end

    % Calculate heated volume
    heated_volume = sum(m_dot_he / rho_w);

    disp(['Hot water produced: ' num2str(heated_volume)])
    if heated_volume >= volume_min
        increase_area = false;
    end
end
disp(Ac)
