% S-calculation

% Year,Month,Day,Hour,Minute,DHI,Temperature,Clearsky DHI,Clearsky DNI,Clearsky GHI,Cloud Type,Dew Point,DNI,GHI,Relative Humidity,Solar Zenith Angle,Surface Albedo,Pressure,Precipitable Water,Wind Direction,Wind Speed
%   1    2    3   4     5     6       7           8             9           10          11        12     13   14        15                 16              17           18           19               20            21 
data = extract_data();
DHI = data(:,6);% Gd
DNI = data(:,13);% Gb
rho_g = data(:,17);%ground albedo
zenith_sun = data(:,16);%solar zenith angle theta_z
beta = data(12,16);%collector zenith angle at 12 noon (max DHI)

n_day = day(datetime(2021,6,21),'dayofyear');%day of the year

% %Given values
[n_c,kl,n_ri,L_back,k_back,L_cover,W,len_tube,diam_tube,len_collector,C_b,L_plate,U,eps_req] = given_vals();

% 47.6205° N, 122.3493° W
latitude = 47.6205;
longitude = 122.3493;
azimuth = 0;% collector azimuth straight south
Gsc = 1367; % W/m^2


declination = 23.45*sind(360*(284+n_day)/365);% 1.6.1a, delta
% omega = [-15*11:15:15*11]';% hour angle (15deg per hour off solar noon)
% omega1_rise = omega(4);
% omega2_rise = omega(5);
% omega1_set = omega(19);
% omega2_set = omega(20);

% cosd_incidence_angle2 = sind(declination).*sind(latitude).*cosd(beta)...
%                      - sind(declination).*cosd(latitude).*sind(beta).*cosd(azimuth)...
%                      + cosd(declination).*cosd(latitude).*cosd(beta).*cosd(omega)...
%                      + cosd(declination).*sind(latitude).*sind(beta).*cosd(azimuth).*cosd(omega)...
%                      + cosd(declination).*sind(beta).*sind(azimuth).*sind(omega);%1.6.2

azimuth_s = sign(omega).*abs(acosd((cosd(zenith_sun).*sind(latitude)-sind(declination))...
                                        ./(sind(zenith_sun).*cosd(latitude)))); %1.6.6

cosd_incidence_angle3 = cosd(zenith_sun)*cosd(beta) + sind(zenith_sun).*sind(beta).*cosd(azimuth_s - azimuth);%1.6.3 TODO: compare with 1.6.2
incidence_angle = acosd(cosd_incidence_angle3);
cosd_zenith_sun = cosd(latitude).*cosd(declination).*cosd(omega) + sind(latitude).*sind(declination);%1.6.5
Go = Gsc*(1 + 0.033*cosd(360*n_day/365).*cosd(zenith_sun));%1.10.1

Id = DHI*3600;%diffuse
Ib = DNI*3600;%beam
Io = Go*3600;%extraterrestrial
I = Ib + Id; %TODO: check

Ai = Ib./Io;% 2.16.3 anisotropy index
f = sqrt(Ib./I);%2.16.6

%% Ratio between incident radiation on horizontal and tilted surface
% R_b1 = cosd_incidence_angle3./cosd(zenith_sun); %1.8.1
R_b2 = (cosd(latitude-beta).*cosd(declination).*cosd(omega) + sind(latitude-beta).*sind(declination))...
        ./ cosd(latitude).*cosd(declination).*cosd(omega) + sind(latitude).*sind(declination);% 1.8.2 TODO: Doesn't match with 1.8.1

% Account for sun rise/set
% a_rise = (sind(declination).*sind(latitude).*cosd(beta) - sind(declination).*cosd(latitude).*sind(beta).*cosd(azimuth)).*(omega2_rise-omega1_rise).*pi/180 ...
%        + (cosd(declination).*cosd(latitude).*cosd(beta) + cosd(declination).*sind(latitude).*sind(beta).*cosd(azimuth)).*(sind(omega2_rise)-sind(omega1_rise)) ...
%        - (cosd(declination).*sind(beta).*sind(azimuth)).*(cosd(omega2_rise)-cosd(omega1_rise));
% b_rise = (cosd(latitude).*cosd(declination)).*(sind(omega2_rise)-sind(omega1_rise)) ...
%         + sind(latitude).*sind(declination).*(omega2_rise-omega1_rise).*pi/180;
% Rb_rise = a_rise/b_rise; %2.14.6
% 
% a_set = (sind(declination).*sind(latitude).*cosd(beta) - sind(declination).*cosd(latitude).*sind(beta).*cosd(azimuth)).*(omega2_set-omega1_set).*pi/180 ...
%        + (cosd(declination).*cosd(latitude).*cosd(beta) + cosd(declination).*sind(latitude).*sind(beta).*cosd(azimuth)).*(sind(omega2_set)-sind(omega1_set)) ...
%        - (cosd(declination).*sind(beta).*sind(azimuth)).*(cosd(omega2_set)-cosd(omega1_set));
% b_set = (cosd(latitude).*cosd(declination)).*(sind(omega2_set)-sind(omega1_set)) ...
%         + sind(latitude).*sind(declination).*(omega2_set-omega1_set).*pi/180; %2.14.6
% Rb_set = a_set./b_set;

Rb = R_b2;
% incorporate sunset corrections TODO: Corrections are very low
% Rb(5) = Rb_rise;
% Rb(19) = Rb_set;



IT = (Ib+Id.*Ai).*Rb + Id.*(1-Ai).*(1+cosd(beta))*0.5.*(1+f.*sind(beta/2).^3) + I.*rho_g.*(1-cosd(beta))*0.5;%2.16.7

%% Absorbtance and transmittance calcs
%from HW2 P1, assume T = 400K 
alpha_n = 0.944;%fig 4.8.3 plate absorbtance
eps = 0.095;%4.8.3 plate emissivity

%fig. 5.4.1
theta_ed = 58;%effective angle of incidence for diffuse sky radiation
theta_eg = 76.5;%ground reflected radiation

tau_alpha_n = 1.01*0.83*alpha_n;%5.5.2, fig 5.3.1

tau_b = [0;0;0;0;0;0.17;0.55;0.75;0.82;0.83;0.83;0.83;0.83;0.83;0.82;0.77;0.61;0.28;0;0;0;0;0];%fig 5.3.1 transmittance from solar incidence angle

%fig 5.6.1
tau_alpha_d = tau_alpha_n*0.84;
tau_alpha_g = tau_alpha_n*0.36;
tau_alpha_b = tau_alpha_n*[0;0;0;0;0;0.1;0.55;0.82;0.95;0.98;0.99;1;0.99;0.98;0.96;0.89;0.64;0.18;0;0;0;0;0];

S = Ib.*Rb.*tau_alpha_b + Id.*tau_alpha_d.*(1+cosd(beta))*0.5 + rho_g.*I.*tau_alpha_g.*(1-cosd(beta)).*0.5;%5.9.1




