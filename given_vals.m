%Run at the top of every script to import constants

function [n_c,kl,n_ri,L_back,k_back,L_cover,W,len_tube,diam_tube,len_collector,C_b,L_plate,U,eps_req] = given_vals()
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
    C_b = 1e6;%bond conductance between tubes and absorber
    L_plate = 5e-4;%absorber plate thickness
    U = 500;%overall heat exchanger coefficient
    eps_req = 0.85;%required effectiveness at every hour

    %water properties found here:

    %copper properties found here:

end