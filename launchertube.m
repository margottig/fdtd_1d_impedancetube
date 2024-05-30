clear all;
close all;
clc;

dh=.01; %spatial step size (1 cm)
c = 341; 
dt = dh/max(c);
iimax=15/341/dt; % simulation will run for iimax iterations


%% Given values
rho1 = 1.21;   % Density of the first medium (kg/m^3)
c1 = 341;      % Speed of sound in the first medium (m/s)
rho2 = 40*1.21;   % Density of the second medium (kg/m^3)
c2 = 341;      % Speed of sound in the second medium (m/s)
L = 1;       % Length of the second medium (meters)
z1 = rho1 * c1;
z2 = rho2 * c2;

%% CALL THE FUNCTION - Record impulse response
%store the sound pressure at two different points in the impedance tube over time
[h1, h2]=first_assignment(1,0.01/341);

%% Get the FFT
%Explain why nonzeros and also plot the whole spectrum
fft_h1 = abs(fft(h1));
fft_h2 = abs(fft(h2));

%% Ricker Wavelet - EXITATION OF THE TUBE
% figure(1);
% plot(h1)
% hold on
% plot(h2);
% title("Magnitude of Ricker Wavelet");
% xlabel('steps');
% ylabel('Magnitude');
% legend('Excitation  signal', 'Transmitted signal');

%% Ricker Wavelet - TIME DOMAIN
figure(2);
plot((1:iimax)*.01/341, h1)
hold on
plot((1:iimax)*.01/341, h2);
title("Magnitude of Ricker Wavelet");
xlabel('Time Domain [s]');
ylabel('Magnitude');
legend('Excitation  signal', 'Transmitted signal');


%% PLOT THE SPECTRUM of the WAVELET in the frequency domain
figure(3)
n = length(h1); f=[0:n/2-1 -n/2:-1]/n*(341/0.01);
plot(f,abs(fft(h1)), f,abs(fft(h2)))
title("Spectrum of the Ricker wavelet vs. frequency (Hz)");
xlabel('Frequency (Hz)');
ylabel('Sound pressure');
legend('Direct Sound', 'Transmitted sound');

%% Energy Transmission Coefficient
% WINDOWING
h = h1;
h(500:end) = 0;
n= 2^14;
dt2 = 0.01/341;
freq = (0:n/2-1)/n/dt;

H1 = abs(fft(h,n)); H1=H1(1:end/2);
H2 = abs(fft(h2,n)); H2=H2(1:end/2);

% Theoretical Value
theoretical = 4./(4+(((z2/z1)-(z1/z2))^2).*((sin(2*pi*freq/341)).^2));

figure(5)
plot(freq, H1.^2, freq, H2.^2);
title("Incident Energy vs. Transmitted Energy");
xlabel('Frequency (Hz)');
ylabel('Sound pressure');
legend('Direct Sound', 'Transmitted sound');

figure(6)
plot(freq,H2.^2./H1.^2, freq, theoretical );
xlim([0 5000])
title("Energy Transmission Coefficient");
xlabel('Frequency (Hz)');
ylabel('Sound pressure');
legend('Measured', 'Theoretical');



