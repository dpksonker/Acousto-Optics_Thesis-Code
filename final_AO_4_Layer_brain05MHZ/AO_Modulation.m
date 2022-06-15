function [Photon,Tissue,M0,M1,AO,mod_dem,phi_dem] =  AO_Modulation(Photon,Tissue,Pn_cosi,Pn_sini,Pd_cosj,Pd_sinj,Bound,M0,M1,dwa,AO,mod_dem,phi_dem)
%% The function AO_Modulation  calculates the total modulation magnitude and phase of the phase modulation and also stores the magnitude 
%% and phase changed in each layers
 
% INPUTS:

% Photon:  Structure to photon variables.
% Bound:   Structure to medium dimensions.
% Tissue:  Structure to tissue variables.
% Pn_cosi: Temporary varible for storing the phase accumulations due to
%          change in refractive index
% Pn_sini: Temporary varible for storing the phase accumulations due to
%          change in refractive index
% Pn_cosj: Temporary varible for storing the phase accumulations due to
%          scatterer displacement
% Pn_sinj: Temporary varible for storing the phase accumulations due to
%          scatterer displacement
% M0:      Variable to accumulate the zeroth order modulation
% M1:      Variable to accumulate the first order modulation
% dwa :    The dropped weight of the photon
% AO :     Structure to loaded variable matrix of Pressure magnitude, Phase
%         ,and direction of US
% mod_dem : Variable for storing the magnitude of phase modulation in each
%           layers
% phi_dem : Variable for storing the phase of phase modulation in each
%           layers


% OUTPUTS:

% Tissue: Structure to tissue variables.
% Pn_cosi: Temporary varible for storing the phase accumulations
% Pn_sini: Temporary varible for storing the phase accumulations
% Photon: Structure containing all the variables that photon carries
% AO :    Structure to loaded variable matrix of Pressure magnitude, Phase
%         , and direction of US
% M0:      Variable to accumulate the zeroth order modulation
% M1:      Variable to accumulate the first order modulation
% mod_dem : Variable for storing the magnitude of phase modulation in each
%           layers
% phi_dem : Variable for storing the phase of phase modulation in each
%           layers

% % The applied equation is from the article: [Sakadžić, Sava, and Lihong
% V. Wang. "Correlation transfer equation for ultrasound-modulated multiply scattered light." Physical Review E 74.3 (2006): 036618.]

%%

[Layer] = Layer_P(Photon,Tissue);                                                % Retrieve the current layer of the photon  

ixd = round(abs(Photon.x/Bound.dx));                                             % Calculate the index of photon in x

if ixd == 0
    
    ix = single(1);
    
elseif ixd > (Bound.Nx-1)

    ix = Bound.Nx-1;                                                                % if not then take the last index   

else
    ix = ixd;                                                                       % if yes then take the converted index 
    
end

iyd = round(abs(Photon.y/Bound.dy));                                                 % Calculate the index of photon in y

if iyd == 0
    
    iy = single(1);
    
elseif iyd > (Bound.Ny-1)
    
    iy = Bound.Ny-1;
    
else
    iy = iyd;
    
end

% 
izd = round(abs(Photon.z/Bound.dz));                                                 % Calculate the index of photon in z

if izd == single(0)

    iz = single(1);
elseif izd > (Bound.Nz-1)

    iz = Bound.Nz-1;

else
    iz = izd;

end


mod_phi = sqrt(Pn_cosi*Pn_cosi + Pd_cosj*Pd_cosj + Pn_sini*Pn_sini + Pd_sinj*Pd_sinj);      % Calculate the magnitude of the modulation
phi     =  (- Pn_sini - Pd_sinj)/(Pn_cosi + Pd_cosj);                                       % Calculate the phase of the modulation

mod_dem(1,Photon.Layer) = mod_dem(1,Photon.Layer) + mod_phi;                                % Store the modulation magnitude in a specific layer
phi_dem(1,Photon.Layer) = phi_dem(1,Photon.Layer) + phi;                                    % Store the modulation phase in a specific layer
% 
% 
M0(ix,iy,iz) = M0(ix,iy,iz) + real(besselj(0,abs(mod_phi))*besselj(0,abs(mod_phi)))*dwa/Layer.mu_a;     % Calculate the zeroth order modulation

M1(ix,iy,iz) = M1(ix,iy,iz) + 2*real(besselj(1,abs(mod_phi))*besselj(1,abs(mod_phi)))*dwa/Layer.mu_a;    % Calculate the first order modulation



end
