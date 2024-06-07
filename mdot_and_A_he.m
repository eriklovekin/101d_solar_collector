function mdot_he,maxA = mdot_and_A_he(eps,mdot_sc);

    mdot_he=eps.*(61.3-20)./35; % from 3.17.5
    U=500; %given in problem statement
    Rc=mdot_sc./mdot_he;
    Ntu = 1./(1-Rc).*ln((1-eps.*Rc)./1-eps); % from Table 8.3b
    A=Ntu.*mdot_he.*4180./U; % from 3.17.8; 4180 from lookup table
    
    return mdot_he,max(A)

end