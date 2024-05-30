function [h1, h2]=first_assignment(pos_1, dt)
% Calculate the impulse response in a specific point in an impedance tube for a given time step.
% [h1 h2] = first_assignmentMA(dt, pos_1)
% dt Temporal step size (in seconds) ~ .01/341;

%% Spatial Discretization
dh=.01; %spatial step size (1 cm)
nx=round(10/dh); % nx= how many steps,  size of the domain = nx*dh (10 m)


%% Simulation Parameters
rho1=1.21;c1=341; k1=(c1^2)*rho1; %first media constants
rho2=48.4;c2=341;k2=(c2^2)*rho2; %second media constants
rho3=1.21; c3=341; k3=(c3^2)*rho3; % Third media constants
dt = dh/max(c1,c2);
iimax=15/341/dt; %simulation will run for iimax iterations

%% Sizes calculation for each medium segment ensuring they sum to nx
% Since dividing nx by 3 might not yield an integer, especially after rounding, calculate integer values for the segment sizes of each medium.
n1 = floor(nx * 0.5);  % 50% for the first medium
n2 = floor(nx * 0.1);  % 10% for the second medium (absorptive material)
n3 = nx - n1 - n2;  % 40% for the third medium, adjusted to fit the remaining length

%% Constructing k and rho arrays
k = [repmat(k1, 1, n1), repmat(k2, 1, n2), repmat(k3, 1, n3)].';
rho = [repmat(rho1, 1, n1-1), (rho2+rho1)/2, repmat(rho2, 1, n2-1), (rho3+rho2)/2, repmat(rho3, 1, n3-1)].';

%% Initialization
p=zeros(nx,1);    %sound preasure
ux=zeros(nx+1,1); %particle velocity
h1=[];      % source impulse response 
h2=[];      % load impulse response

%% Exitation Setup (vectorized operations) 
%excitation
centralfrequency=1000;
a=8*centralfrequency/sqrt(pi);
t=((1:iimax)/(1/dt)-4/a);  %time vector  t that is used to evaluate an excitation signal w.
w=-(exp(-a^2*(t.^2)/2).*(a^2*(t.^2)-1)); %Ricker wavelet (mexican hat)

%% FDTD operation
for ii=1:iimax
    p=p-k.*dt/dh.*(diff(ux));                  %preassure calculation
    ux(2:nx,1)=ux(2:nx,1)-dt./rho./dh.*diff(p); %velocity calculation
  
    %% Boundary Conditions
    ux(end,1)=p(end)/rho3/c3;       %boundary conditions (rho3, so impedance = 1, then total absorption) if using rho2 we have reflections instead4444444
    ux(1)=w(ii)/c1/rho1;    %excitation (beggining 1st medium)
   
    h1=[h1 p(round(pos_1))]; % Record the pressure near the sound source of the entire domain 
    h2=[h2 p(round(601))]; % Record the pressure once it had pass the load (absorpive material)  

    %% PLOT IN TERMS OF SOUND PRESSURE 
    % blue_line= pressure, red_line= particle velocity
    %plot((1:nx)*dh,10*log10(max(p.^2/(2e-5).^2,1)),(.5:nx+.5)*dh,10*log10(max(ux.^2/(2e-5/413).^2,1))) 
    %axis([0 10 20 100]);drawnow
end
return