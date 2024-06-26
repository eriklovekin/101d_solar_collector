%Run at the top of every script to import constants

function [n_c,kl,n_ri,L_back,k_back,L_cover,W,len_tube,diam_tube,len_collector,C_b,L_plate,U,eps_req,cp_w,rho_w,mu_w,nu_w,Pr_w,k_w,k_c] = given_vals()
    % Given in problem statement
    n_c = 2;% number of covers
    kl = 0.0125;% cover glass KL product (per sheet)
    n_ri = 1.526;% cover refractive index
    L_back = 0.007;%back insulation thickness
    k_back = 0.0245;%back insulation thermal conductivity
    L_cover = 0.025;%cover spacing
    W = 0.115;%distance between tubes
    len_tube = 2.46;%tube length
    diam_tube = 0.015;%tube diameter
    len_collector = 2.5;%collector gross length
    C_b = 1e7;%bond conductance between tubes and absorber
    L_plate = 5e-4;%absorber plate thickness
    U = 500;%overall heat exchanger coefficient
    eps_req = 0.85;%required effectiveness at every hour

    %material properties of water at 300K from Basic Heat Transfer Table A.8
    cp_w = 4178;% J/kgK
    rho_w = 996;%kg/m^3
    mu_w = 8.67e-4;%kg/m/s
    nu_w = 0.87e-6;%m^2/s
    Pr_w = 5.9;
    k_w = 0.611;%W/mK

    %copper properties found here: https://www.matweb.com/search/datasheet.aspx?matguid=9aebe83845c04c1db5126fada6f76f7e&ckck=1
    k_c = 398;% W/mK

end