clc
clearvars


load('Data_05.mat')

% Mag_mod1 = Data_05MHzFUS_4LBRAIN{1, 1}(:,2:end);
Mag_mod1 = Data_05MHzFUS_water{1, 1}(:,2);
Mag_mod2 = Data_05MHzFUS_water{2, 1}(:,2);
Mag_mod3 = Data_05MHzFUS_water{3, 1}(:,2);

phi_mod1 = Data_water_05MHz{4, 1}(:,2);
phi_mod2 = Data_water_05MHz{5, 1}(:,2);
phi_mod3 = Data_water_05MHz{6, 1}(:,2);

ppl1 = Data_water_05MHz{7, 1}(:,3); 
ppl2 = Data_water_05MHz{8, 1}(:,3); 
ppl3 = Data_water_05MHz{9, 1}(:,3); 

mu_a = 0.035;
dt = (1e-3/8.5)/(3e8);

EnergyD1 = exp(-mu_a*ppl1);

Flux0 = EnergyD1.*(besselj(0,Mag_mod1).^2);
Flux1 = 2*EnergyD1.*(besselj(1,Mag_mod1).^2);

MD = Flux1./Flux0;
figure(1)

plot(Flux0)
hold on
plot(Flux1,'r')
hold on
plot(MD)

figure(2)
plot(Flux0/max(Flux0))
hold on
plot(Flux1/max(Flux1))
hold on
plot(MD/max(MD))


n = length(Flux0);
fs = 1/dt;
XX = fftshift(fft(abs(MD)));
freq1 = (-n/2:n/2-1)*(fs/n);

power = abs(XX).^2/n;

figure(34)
plot(freq1, power)
% 