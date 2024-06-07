%

Tci = 20;
Tco = 55;

Thi = 61.2;

volume_min = 120*0.00378541;%m^3
% Year,Month,Day,Hour,Minute,DHI,Temperature,Clearsky DHI,Clearsky DNI,Clearsky GHI,Cloud Type,Dew Point,DNI,GHI,Relative Humidity,Solar Zenith Angle,Surface Albedo,Pressure,Precipitable Water,Wind Direction,Wind Speed
%   1    2    3   4     5     6       7           8             9           10          11        12     13   14        15                 16              17           18           19               20            21 
[n_c,kl,n_ri,L_back,k_back,L_cover,W,len_tube,D,len_collector,Cb,delta,U,eps_req,cp_w,rho_w,mu_w,nu_w,Pr_w,k_w,k_c] = given_vals();
data = extract_data();
wind_speed = data(:,21);
wind_dir = data(:,20);
Ta = data(:,7);


n = 0;% number of tubes
increase_area = true;
conv = 0.01;

while increase_area
    n = n+1;
    Ac = len_collector*n*len_tube;
    
       
    
        % Assuming Ch = Cmin, 3.17.5 simplifies
        Tho = Tco - eps_req*(Thi-Tci);

        %Guess values to calculate m_dot_h
        m_dot_h = 0.1;%kg/s
        Tpm = (Thi+Tho)/2;
        Tc1 = Tpm;
        Tc2 = Tpm;

        UL = UL_calc(Tc1,Tc2,Tpm);

        hfi = hfi_calc(m_dot,Thi,Tho);% convective heat transfer coefficient inside tube (doesn't use Thi,Tho currently)

        m = sqrt(UL/(k_c*delta));% DB 6.5.4a
        F = tanh((m*(W-D)/2)/(m*(W-D)/2));% DB 
        
        
        tmp1 = 1/(UL*(D + (W - D)*F));
        tmp2 = 1/(pi*D*hfi);
        F_prime = (1/UL)/(W*(tmp1 + 1/Cb + tmp2));










    if heated_volume >= volume_min
        increase_area = false;
    end
end




