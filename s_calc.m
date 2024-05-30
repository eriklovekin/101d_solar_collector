% S-calculation

% Year,Month,Day,Hour,Minute,DHI,Temperature,Clearsky DHI,Clearsky DNI,Clearsky GHI,Cloud Type,Dew Point,DNI,GHI,Relative Humidity,Solar Zenith Angle,Surface Albedo,Pressure,Precipitable Water,Wind Direction,Wind Speed
%   1    2    3   4     5     6       7           8             9           10          11        12     13   14        15                 16              17           18           19               20            21 
data = extract_data();
DHI = data(:,6);%
T_a = data(:,7);%air temp
theta_z = data(:,16);%solar zenith angle
rho_g = data(:,17);%ground albedo
dir_wind = data(:,20);%wind direction
v_wind = data(:,21);%wind speed

n_day = day(datetime(2021,6,21),'dayofyear');%day of the year

%Given values
[n_c,kl,n_ri,L_back,k_back,L_cover,W,len_tube,diam_tube,len_collector,C_b,L_plate,U,eps_req] = given_vals();

R_b = 


% Anisotropy index
A_i = I_bn/I_on = I_b/I_o;% 2.16.3 
f = sqrt(I_b/I);% 2.16.6

% Total radiation on tilted surface using Reindl et al.
I_dT = I_d*((1-A_i)*((1 + cosd(beta))/2)*(1 + f*sind(beta/2)^3) + A_i*R_b);%2.16.5
I_T = (I_b + I_d*A_i)*R_b + I_d*(1 - A_i)*((1+cosd(beta))/2)*(1 + f*sind(beta/2)^3) + I*rho_g*((1-cosd(beta))/2);% 2.16.7



