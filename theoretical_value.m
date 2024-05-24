function [theo_transmission_coefficients]=theoretical_value(numberOfPoints)
% Given values
rho1 = 1.21;   % Density of the first medium (kg/m^3)
c1 = 341;      % Speed of sound in the first medium (m/s)
rho2 = 48.4;   % Density of the second medium (kg/m^3)
c2 = 341;      % Speed of sound in the second medium (m/s)
L = 1;       % Length of the second medium (meters)

% Frequency range (from 0 to 300 Hz)
frequencies = linspace(0, 34100/2 ,numberOfPoints); %Oversampling

% Initialize an array to store transmission coefficients
theo_transmission_coefficients = zeros(size(frequencies));

% Calculate acoustic impedances
z1 = rho1 * c1;
z2 = rho2 * c2;

% Loop over each frequency and calculate the transmission coefficient
for ii = 1:length(frequencies)
    f = frequencies(ii);
    
    % Calculate wave number in the second medium
    k2 = 2 * pi * f / c2;
    
    % Calculate the sine term
    sin_term = sin(k2 * L).^2;
    
    % Calculate the transmission coefficient
    numerator = 4;
    denominator = 4 + ((((z2 / z1) - (z1 / z2)).^2 ).* sin_term);
   theo_transmission_coefficients(ii) = numerator ./ denominator;
end



