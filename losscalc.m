data = extract_data();

n_covers = 2;
KL_covers = .0125;
RI_covers = 1.526;
% wind speed coefficcient needs to be calced
windspeeds = data(:,end); % but these will vary based on angle of collector (in m/s presumably)
back_insulation_thickness = .007;
back_insulation_conductivity = .0245;
spacing = 0.025; 
tube_distance = 0.115;
tube_length = 2.46;
tube_diameter = 0.015; % thin walled
collector_length = 2.5;
bond_conductance = 10e6;
plate_thickness = 0.0005;
