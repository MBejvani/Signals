clc , clear
%% SignalClass:
% A class encompassing various types of signals
% Written by Mohammad Bejvani
% Email--> mohammadbejvani@gmail.com and in Github--> @MBejvani
%% S Groups
fs = 100; dt = 1/fs;
t = 0:dt:10-dt;
N = length(t);
f = linspace(0,1/dt/2,N/2);
Signals.S.dt = dt;
% MonoFreq3C
fi_I = 2/5*(t-t(round(N/2))).^3 + 3*(t-t(round(N/2))) ;
s1 = sin(2*pi*f(round(N/10))*t) + sin(2*pi*f(round(N/4))*t) + sin(2*pi*f(round(3*N/8))*t);
Signals.S.MonoFreq3C = s1/max(s1);
% FullFreq3C
Signals.S.FullFreq3C = zeros(1,N);s2(round(N/5)) = 6;s2(round(N/2)) = 6;s2(round(4*N/5)) = 6;
% HypChirps
Signals.S.HypChirps = sin(2*pi*fi_I);
% Wavelets
Signals.S.Wavelets = sin(2*pi*f(round(9*N/20))*t.*(.8<t & t<1)) + ...
    sin(2*pi*f(round(9*N/20))*t.*(3.3<t & t<3.5)) + ...
    sin(2*pi*f(round(9*N/20))*t.*(6.4<t & t<6.6)) + ...
    sin(2*pi*f(round(9*N/20))*t.*(8.9<t & t<9.1));
% TwoWavelets
Signals.S.TwoWavelets = sin(2*pi*f(round(N/6))*t) .* ...
    exp(-1/2*(((0:N-1)-fix(N/3))/35).^2) + sin(2*pi*f(round(N/6))*t) .* ...
    exp(-1/2*(((0:N-1)-fix(2*N/3))/35).^2);
%% Sinusoidal-IF Siganl 
N = 1024;
t = linspace(0,1,N);
phi1 = 5 * sin(15*t)+ 200 * t;
cr = 5*15*-15*sin(15*t);
s = cos(2*pi*phi1);
dt = t(2)-t(1);
Signals.SinIF.dt=dt;
Signals.SinIF.s=s;
%% 2-component Sinusoidal-IF signal  
N = 1024;
t = linspace(0,1,N);
phi1 = 40*pi*t-10*t.^2; % Phase
phi2 = 5*sin(25*t)+100*pi*t; % Phase
IF0 = (125*cos(25*t)+130*pi-20*t)/2; % Instantaneous Frequency
CR0 = (-125*25*sin(25*t)-20)/2; % Chirp-rate
s = .5*cos(2*pi*phi1) + .5*cos(2*pi*phi2) ;
dt = t(2)-t(1);
f = linspace(-1/dt/2,1/dt/2,N);
Signals.SinIF2C.dt=dt;
Signals.SinIF2C.s=s;
Signals.SinIF2C.t=t;
Signals.SinIF2C.f=f;
Signals.SinIF2C.cr=CR0;
%% Linear Chip
N = 1024;
t = linspace(0,1,N);
dt = t(2)-t(1);
phi1 = 1/(4*dt*max(t))*t.^2;
s = cos(2*pi*phi1);
Signals.LinChirp.dt=dt;
Signals.LinChirp.s=s;
%% Three Stairs 
N = 1024;
t = linspace(0,1,N);
dt = t(2)-t(1);
s = [cos(2*pi*t(1:fix(N/3))*100) cos(2*pi*t(1+fix(N/3):fix(2*N/3))*200) cos(2*pi*t(1+fix(2*N/3):end)*300)];
Signals.Stairs.dt=dt;
Signals.Stairs.s=s;
%% Monofrequency
N = 1024;
t = linspace(0,1,N);
dt = t(2)-t(1);
phi1 = 1 / (5 * dt ) * t;
s = cos(2*pi*phi1);
Signals.MonoFreq.dt=dt;
Signals.MonoFreq.s=s;
%% 3-component exponential chirps
N = 1024;
t = linspace(0,1,N)+eps;
phi1 = 1/6*exp(6*t)+25*t;
phi2 = 1/6.2*exp(6.2*t)+50*t;
phi3 = 1/5.8*exp(5.8*t)+t;
cr = 6*exp(6*t);
s = cos(2*pi*phi1) + cos(2*pi*phi2) + cos(2*pi*phi3);
dt = t(2)-t(1);
Signals.ExpChip3C.dt=dt;
Signals.ExpChip3C.s=s;
Signals.ExpChip3C.cr=cr;
%% 2-component exponential chirps
N = 1024;
t = linspace(0,1,N)+eps;
phi1 = 1/6*exp(6*t)+25*t;
phi2 = 1/6.2*exp(6.2*t)+50*t;
s = cos(2*pi*phi1) + cos(2*pi*phi2);
dt = t(2)-t(1);
Signals.ExpChip2C.dt=dt;
Signals.ExpChip2C.s=s;
%% Droped
N = 1024;
t = linspace(0,1,N)+eps;
phi1 = -300*(t(fix(N/2)+1:end)-t(fix(N/2))).^3 + 400 * t(fix(N/2)+1:end);
phi2 = 300*(t(1:fix(N/2))-t(fix(N/2))).^3 + 100 * t(1:fix(N/2));
cr1 = 2*3*-300*(t(fix(N/2)+1:end)-t(fix(N/2)));
cr2 = 2*3*300*(t(1:fix(N/2))-t(fix(N/2)));
cr = [cr1 cr2];
s = [cos(2*pi*phi1) cos(2*pi*phi2)];
dt = t(2)-t(1);
Signals.Droped.dt=dt;
Signals.Droped.s=s;
%% RickerWavelet
dt = .008 ;
trng = 8 ; 
df = 1/trng ;
t = 0:dt:trng-dt ;
N = fix(trng/dt);
r([100,300,350,490,500,650,660,670,860,880,900]) = [2 1.5 -2 1 1 -.75 .75 -.75 .5 -.5 .5];
fd = 10;
x = (1-2*(pi*fd*(t-trng/2)).^2) .* exp(-(pi*fd*(t-trng/2)).^2) ;
y = conv(r,x);s = y(fix(N/2)+1:N+fix(N/2));
s = s + .5*cos(2*pi*25*t);
Signals.RickerWavelet.dt=dt;
Signals.RickerWavelet.s=s;
%% Trace Signal
load('trace.mat')
Signals.RealTrace=trace;
%% Save Signals
save ('SignalsTypes','Signals')