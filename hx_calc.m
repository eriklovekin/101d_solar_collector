%

Tci = 20;
Tco = 55;

Thi = 61.2;

volume_min = 120*0.00378541;%m^3

[n_c,kl,n_ri,L_back,k_back,L_cover,W,len_tube,diam_tube,len_collector,C_b,L_plate,U,eps_req,cp_w,rho_w,mu_w,nu_w,Pr_w] = given_vals();
data = extract_data();
Ta = data(:,7);


n = 0;% number of tubes
increase_area = true;
conv = 0.01;

while increase_area
    n = n+1;
    Ac = len_collector*n*len_tube;
    for hour = 5:19
        m_dot_h = 0.1;%kg/s initial guess

        % Assuming Ch = Cmin, 3.17.5 simplifies
        Tho = Tco - eps_req*(Thi-Tci);

        %Guess values to calculate m_dot_h
        m_dot_h = 0.1;%kg/s
        Tpm = (Thi+Tho)/2;
        Tc1 = Tpm;
        Tc2 = Tpm;
        while abs(m_dot_h - m_dot_h_old) >= conv
            %green while loop

            
            Re = 4*m_dot_h/(pi*diam_tube*mu_w);%3.14.1
            
            hfi = 

        end






    end
    if heated_volume >= volume_min
        increase_area = false;
    end
end




