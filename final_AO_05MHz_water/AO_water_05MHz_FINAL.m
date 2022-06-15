clearvars
clc

tic
% Tissue Parameters
load('Pressure.mat');
load('Phase.mat');
load('u_x.mat');
load('u_y.mat');
load('u_z.mat');

AO.Pressure = n_p/5.3837e3;                                                 % normalize the US pressure to 1 kPa
AO.Phase = n_ph;                                                            % upload the phase data
AO.ux = u_xx;                                                               % particle velocity in x direction
AO.uy = u_yy;
AO.uz = u_zz;

Bound.Nx = single(202);
Bound.Ny = single(202);
Bound.Nz = single(101);                                                            % Dimension of the tissue


Bound.dx = single(0.5);
Bound.dy = single(0.5);                                                           % Spacing between grid lines 0.5 mm
Bound.dz = single(0.5);


Bound.x = single(zeros(1,202));
Bound.y = single(zeros(1,202));                                                  % Initialize the x,y,z,r to zeros
Bound.z = single(zeros(1,101));


Bound.xmin = single(0);
Bound.xmax = single(0);                                                            % Variable initialized to and are for storing min and max tissue dimensions
Bound.ymin = single(0);
Bound.ymax = single(0);
Bound.zmin = single(0);
Bound.zmax = single(0);

Bound.x = ([0:Bound.Nx-1] - Bound.Nx/2)*Bound.dx;
Bound.y = ([0:Bound.Ny-1] - Bound.Ny/2)*Bound.dy;                          % The coordinates of x,y,z,r
Bound.z = ([0:Bound.Nz-1])* Bound.dz;


Bound.xmin = min(Bound.x);
Bound.xmax = max(Bound.x);
Bound.ymin = min(Bound.y);                                                 % Min and Max tissue length
Bound.ymax = max(Bound.y);
Bound.zmin = min(Bound.z);
Bound.zmax = max(Bound.z);

Tissue.Radius_det = single(0.15);                                                                  % Radius of the detector

Tissue.D1 = single(0);
Tissue.D1 = single(20);                                            % Distance of D1 from source

Tissue.D2 = 0;
Tissue.D2 = single(30);                                                % Distance of D2 from source

Tissue.D3 = 0;
Tissue.D3 = single(40);                                             % Distance of D3 from source

Tissue.count_D1 = single(0);
Tissue.count_D2 =  single(0);                                                                % Counter to count number of photons reached the detector
Tissue.count_D3 =  single(0);

n = single([1,1.38,1]);                                                         % Row vector representing the refractive indexes of all the layers

Tissue.d1 = single(50.5);                                                             % Tickness of first layer                                                           % Tickness of second layer
lambda = single(850e-9);                                                   % Optical wavelength
Tissue.f0 = single(0.5e6);                                                 % Acoustic Frequency
Tissue.k0 = single(2*pi/lambda);                                           % Optical Wave number

Tissue.eta = single(0.32);                                                 % Acousto-optic coefficient

Tissue.mod1 = single(zeros(1e5,1+(length(n)-2)));                          % Variable for storing magnitude of phase modulation at D1 for each layers
Tissue.mod2 = single(zeros(1e5,1+(length(n)-2)));                          % Variable for storing magnitude of phase modulation at D2 for each layers
Tissue.mod3 = single(zeros(1e5,1+(length(n)-2)));                          % Variable for storing magnitude of phase modulation at D3 for each layers
 
Tissue.phi1 = single(zeros(1e5,1+(length(n)-2)));                          % Variable for storing phase of phase modulation at D1 for each layers
Tissue.phi2 = single(zeros(1e5,1+(length(n)-2)));                          % Variable for storing phase of phase modulation at D2 for each layers
Tissue.phi3 = single(zeros(1e5,1+(length(n)-2)));                          % Variable for storing phase of phase modulation at D3 for each layers

Tissue.Det1 = single(zeros(1e5,2+(length(n)-2)));                          % Variable for storing weight of detected photon, and partial path length in each layers at D1
Tissue.Det2 = single(zeros(1e5,2+(length(n)-2)));                          % Variable for storing weight of detected photon, and partial path length in each layers at D2
Tissue.Det3 = single(zeros(1e5,2+(length(n)-2)));                          % Variable for storing weight of detected photon, and partial path length in each layers at D3
Tissue.Tr = single(0);

M0 = single(zeros(Bound.Nx,Bound.Ny,Bound.Nz));                            % Variable for storing the zeroth order modulation
M1 = single(zeros(Bound.Nx,Bound.Ny,Bound.Nz));                            % Variable for storing the first order modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_photons = 1e5;                                                                         % Total Number of Photons to be simulated
W = single(0.0001);                                                                      % Threshold weight for Roulette

% [Bound,Tissue,AO,M0,M1] = AO_05MHz_water(Bound,Tissue,n,n_photons,W,AO,M0,M1);       % Main calling function to simulate the Monte Carlo
[Bound,Tissue,AO,M0,M1] = AO_05MHz_water_mex(Bound,Tissue,n,n_photons,W,AO,M0,M1);

Data_05MHzFUS_water(1:9,:) = {Tissue.mod1(1:Tissue.count_D1,:),"Tissue.Mod1";Tissue.mod2(1:Tissue.count_D2,:),"Tissue.Mod2";Tissue.mod3(1:Tissue.count_D3,:),"Tissue.Mod3";Tissue.phi1(1:Tissue.count_D1,:),"Tissue.phi1";Tissue.phi2(1:Tissue.count_D2,:),"Tissue.phi2";Tissue.phi3(1:Tissue.count_D3,:),"Tissue.phi3";...
                              Tissue.Det1(1:Tissue.count_D1,:),"Tissue.Det1";Tissue.Det2(1:Tissue.count_D2,:),"Tissue.Det2";Tissue.Det3(1:Tissue.count_D1,:),"Tissue.Det3"};

M0_diagonal =  zeros(Bound.Nx,Bound.Nz);
M1_diagonal =  zeros(Bound.Nx,Bound.Nz);

for ii = 1:Bound.Nz
    M0_diagonal(:,ii) = diag(M0(:,:,ii));
    M1_diagonal(:,ii) = diag(M1(:,:,ii));                                   % Store the diagnal elements of the matrix
end

MD_new = M1_diagonal./M0_diagonal;                                          % Modulation Depth Calculation
%
figure(2)
imagesc(Bound.x,Bound.z,log10(M0_diagonal'));                               % Distribution of zeroth order plot
colormap hot
xlim([-50 0])
ylim([0 49.5])
set(findobj(gcf,'type','axes'),'FontName','Helvetica','FontSize',45,'FontWeight','Bold')
set(gcf,'color','w');
xticks([-50 -40 -30 -20 -10 0]);
xticklabels({'0','10','20','30','40', '50'});
yticks([0 10 20 30 40 49.5]);
yticklabels({'0','10','20','30','40', '50'});
xlabel('x (mm)','FontSize',45)
ylabel('z (mm)','FontSize',45);
colorbar
saveas(gcf,'FUS_Unmodulated.fig')

figure (3)
imagesc(Bound.x,Bound.z,log10(M1_diagonal'));                               % Distribution of first order plot
colormap hot
xlim([-50 0])
ylim([0 49.5])
set(findobj(gcf,'type','axes'),'FontName','Helvetica','FontSize',45,'FontWeight','Bold')
set(gcf,'color','w');
xticks([-50 -40 -30 -20 -10 0]);
xticklabels({'0','10','20','30','40', '50'});
yticks([0 10 20 30 40 49.5]);
yticklabels({'0','10','20','30','40', '50'});
xlabel('x (mm)','FontSize',45)
ylabel('z (mm)','FontSize',45);
colorbar
saveas(gcf,'FUS_Modulated.fig')

figure (4)
imagesc(Bound.x,Bound.z,log10(MD_new'));                                    % Distribution of modulation depth
colormap hot
xlim([-50 0])
ylim([0 49.5])
set(findobj(gcf,'type','axes'),'FontName','Helvetica','FontSize',45,'FontWeight','Bold')
set(gcf,'color','w');
xticks([-50 -40 -30 -20 -10 0]);
xticklabels({'0','10','20','30','40', '50'});
yticks([0 10 20 30 40 49.5]);
yticklabels({'0','10','20','30','40', '50'});
xlabel('x (mm)','FontSize',45)
ylabel('z (mm)','FontSize',45);
colorbar
saveas(gcf,'FUS_MDepth.fig')

ppr = squeeze(n_p(:,51,:))';

figure(5)
imagesc(-flip(Bound.x(1:102)),Bound.z,ppr*1e-6)                             % Plot of cropped US pressure magnitude
xlabel('x (mm)','FontSize',45);
ylabel('z (mm)','FontSize',45);
colormap hot
set(findobj(gcf,'type','axes'),'FontName','Helvetica','FontSize',45,'FontWeight','Bold')
h = colorbar('Ticks',[ 0.1 0.2 0.3 0.4 0.5],'TickLabels',{'0.1','0.2','0.3','0.4','0.5'});
ylabel(h,'P_0 (MPa)','FontSize',45);
saveas(gcf,'CROP_US.fig')

toc

% ppl_ini = Tissue.Det1(1:Tissue.count_D1,3);
% time_arr_p = ppl_ini/3e8;
% 
% t_max = (max(ppl_ini) - min(ppl_ini))/(3e8);
% t_min = min(ppl_ini)/(3e8);
% dt = (t_max - t_min);
% 
% tt = linspace(t_min,t_max,1:Tissue.count_D1);
% % dt =  1.2715e-10;
% fs = 1/dt;
% 
% wt = Tissue.Det1(1:Tissue.count_D1,2);
% phi_mag = Tissue.mod1(1:Tissue.count_D1,2);
% 
% f_light = 3e8/lambda;
% 
% Et = wt.*exp(1i*2*pi*f_light*time_arr_p - phi_mag);
% 
% figure(33)
% plot(abs(Et))
% 
% % ff = fftshift(fft(abs(Et)));
% % figure(34)
% % plot(abs(ff))
% 
% n = length(Et);
% %
% XX = fftshift(fft(abs(Et)));
% freq1 = (-n/2:n/2-1)*(fs/n);
% %
% power = abs(XX).^2/n;
% 
% figure(34)
% plot(freq1, power)
% 
% 
% 
