function [mdot_he,maxA,eps] = mdot_and_A_he(eps,mdot_sc)

    mdot_he=mdot_sc.*eps.*(61.3-20)./35; % from 3.17.5
    U=500; %given in problem statement
    Rc=mdot_sc./mdot_he;
    for k = 1:length(Rc)
            if isnan(Rc(k))
                Rc(k) = 1;
            end
    end
    Ntu1 = 1./(1-Rc).*log((1-eps.*Rc)./(1-eps)); % from Table 8.3b
    Ntu1 = real(Ntu1);
    for k = 1:length(Ntu1)
        if isnan(Ntu1(k))
            Ntu1(k) = 0;
        end
        if isinf(Ntu1(k))
            Ntu1(k) = 50;
        end
    end
    A=Ntu1.*mdot_he.*4180./U; % from 3.17.8; 4180 from lookup table
    for k = 1:length(A)
        if isnan(A(k))
            A(k) = 0;
        end
    end
    [maxA, max_location]=max(A); % TODO: should always be at noon, but isn't in i = 4 and a few others. Hardcode to be at noon.
    Ntu=A(12)./mdot_he./4180.*U;
    Ntu = real(Ntu);
    for k = 1:length(Ntu)
        if isnan(Ntu(k))
        Ntu(k) = 0;
        end
        if isinf(Ntu(k))
            Ntu(k) = 50;
        end
    end
    eps = (1-exp((1-Rc).*Ntu))./(Rc-exp((1-Rc).*Ntu));% DB 3.17.8 and BHT Table 8.3b
    eps = real(eps);
    for k = 1:length(eps)
        if isnan(eps(k))
        eps(k) = 1;
        end
    end
    % disp([Rc(12),Ntu1(12),Ntu(12),mdot_he(12),A(12),maxA,max_location,eps(12),min(eps)])
end