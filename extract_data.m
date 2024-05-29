%Extract Irradiation data from data spreadsheet
% Year,Month,Day,Hour,Minute,DHI,Temperature,Clearsky DHI,Clearsky DNI,Clearsky GHI,Cloud Type,Dew Point,DNI,GHI,Relative Humidity,Solar Zenith Angle,Surface Albedo,Pressure,Precipitable Water,Wind Direction,Wind Speed
%   1    2    3   4     5     6       7           8             9           10          11        12     13   14        15                 16              17           18           19               20            21 
function data = extract_data()
    filename = 'day_data.csv';
    data = readmatrix(filename);
    data = data(3:end,:);
end