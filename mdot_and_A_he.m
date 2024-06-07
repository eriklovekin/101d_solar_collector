function [mdot_he,maxA,eps] = mdot_and_A_he(eps,mdot_sc)

    mdot_he=mdot_sc.*eps.*(61.3-20)./35; % from 3.17.5
    U=500; %given in problem statement
    Rc=mdot_sc./mdot_he;
    Ntu = 1./(1-Rc).*log((1-eps.*Rc)./(1-eps)); % from Table 8.3b
    A=Ntu.*mdot_he.*4180./U; % from 3.17.8; 4180 from lookup table
    maxA=max(A);
    Ntu=maxA./mdot_he./4180.*U;
    eps=(Rc-exp((1-Rc).*Ntu))./(1-exp((1-Rc).*Ntu));% DB 3.17.8 and BHT Table 8.3b
end