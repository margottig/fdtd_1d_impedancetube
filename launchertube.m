clear all;
close all;
clc;

%Dont forget to change iimax to 10000, 

%% Given values
rho1 = 1.21;   % Density of the first medium (kg/m^3)
c1 = 341;      % Speed of sound in the first medium (m/s)
rho2 = 40*1.21;   % Density of the second medium (kg/m^3)
c2 = 341;      % Speed of sound in the second medium (m/s)
L = 1;       % Length of the second medium (meters)


%% CALL THE FUNCTION - Record impulse response
%store the sound pressure at two different points in the impedance tube over time
[h1, h2]=first_assignment(1,0.01/341);

%% Get the FFT
%Explain why nonzeros and also plot the whole spectrum
fft_h1 = abs(fft(h1));
fft_h2 = abs(fft(h2));

%% Ricker Wavelet - EXITATION OF THE TUBE
figure(1);
plot(h1)
hold on
plot(h2);
title("Magnitude of Ricker Wavelet");
xlabel('Time Steps');
ylabel('Magnitude');
legend('Excitation  signal', 'Transmitted signal');

%% PLOT THE SPECTRUM of the WAVELET in the frequency domain
figure(2)
n = length(h1); f=[0:n/2-1 -n/2:-1]/n*(341/0.01);
plot(f,abs(fft(h1)), f,abs(fft(h2)))
title("Spectrum of the Ricker wavelet vs. frequency (Hz)");
xlabel('Frequency (Hz)');
ylabel('Sound pressure (Pa)');
legend('Direct Sound', 'Transmitted sound');

%% Energy Transmission Coefficient
%Data cleaning
p2 = nonzeros(h2); 
p1 = nonzeros(h1); p1 = p1(1:length(p2));
ETC=(abs(fft(p2)).^2./(abs(fft(p1)))).^2;
% figure(6)
% plot(ETC)
ETC=ETC(1:end/2);
frequencies = linspace(0, 34100/2 ,length(ETC)); %Oversampling
[theo_transmission_coefficients]=theoretical_value(length(ETC));

figure(3)
plot(frequencies,ETC/max(ETC) );
hold on
plot(frequencies, theo_transmission_coefficients/max(theo_transmission_coefficients));
title("Energy Transmission Coefficient");
xlabel('Frequency (Hz)');
ylabel('Transmission Coefficient (T)');
legend('Experimental results', 'Theoretical results');
grid on
%hold on
%plot(frequencies, transmission_coefficients)
% Gain=20*log10(abs(fft(h))./abs(fft(h(:,1))));
% Gain=Gain(1:end/2,:);
% figure(1);plot(f,Gain)
%  axis([0 2000 0 40])
% legend(num2str(cases'))


%% GET THE PRESSURE 
% figure(4)
% plot(max(20*log10(abs(fft_h1)),1))
% hold on
% plot(max(20*log10(abs(fft_h2)),1))
% xlabel('Frequency (Hz)');
% ylabel('Sound pressure (Pa)');
% legend('ASDAS', 'BSKDLF');