function [Pd_cosj,Pd_sinj,Photon,Tissue,AO] = Accumulate_ScatterD(Photon,Tissue,Bound,Pd_cosj,Pd_sinj,AO,n)
%% The function Accumulate_ScatterD  accumulates the phase change due to scatterer displacement 
 
% INPUTS:

% Photon: Structure to photon variables.
% Bound:  Structure to medium dimensions.
% Tissue: Structure to tissue variables.
% Pn_cosj: Temporary varible for storing the phase accumulations
% Pn_sinj: Temporary varible for storing the phase accumulations
% AO :    Structure to loaded variable matrix of Pressure magnitude, Phase
%         , and direction of US
% n :     Tissue layers refractive index matrix 

% OUTPUTS:

% Tissue: Structure to tissue variables.
% Pn_cosj: Temporary varible for storing the phase accumulations
% Pn_sinj: Temporary varible for storing the phase accumulations
% Photon: Structure containing all the variables that photon carries
% AO :    Structure to loaded variable matrix of Pressure magnitude, Phase
%         , and direction of US

% % The applied equation is from the article: [Sakadžić, Sava, and Lihong
% V. Wang. "Correlation transfer equation for ultrasound-modulated multiply scattered light." Physical Review E 74.3 (2006): 036618.]
%%
[Layer] = Layer_P(Photon,Tissue);                                               % Retrieve the current layer of the photon         

ixd = round(abs(Photon.x/Bound.dx));                                            % Calculate the index of photon in x

if ixd == 0                                                                     
    
    ix = single(1);
    
elseif ixd > (Bound.Nx-1)                                                       % Check if the converted index is in the bounds of the total index
    
    ix = Bound.Nx-1;                                                            % if not then take the last index   
    
else
    ix = ixd;                                                                   % if yes then take the converted index   
    
end

iyd = round(abs(Photon.y/Bound.dy));                                            % Calculate the index of photon in y

if iyd == 0
    
    iy = single(1);
    
elseif iyd > (Bound.Ny-1)
    
    iy = Bound.Ny-1;
    
else
    iy = iyd;
    
end


izd = round(abs(Photon.z/Bound.dz));                                             % Calculate the index of photon in z

if izd == 0
    
    iz = single(1);
elseif izd > (Bound.Nz-1)
    
    iz = Bound.Nz-1;
    
else
    iz = izd;
    
end
n_layer = Photon.Layer;                                                         % number of layer
Pressure_mag = sqrt(AO.ux(ix,iy,iz)*AO.ux(ix,iy,iz) + AO.uy(ix,iy,iz)*AO.uy(ix,iy,iz) + AO.uz(ix,iy,iz)*AO.uz(ix,iy,iz));      % Calculate the magnitude of US direction vectors

xdiff = Photon.ux_old - Photon.ux;                                                                                             % Calculate the dorection
ydiff = Photon.uy_old - Photon.uy;
zdiff = Photon.uz_old - Photon.uz;

dot_a_o = (AO.ux(ix,iy,iz)/Pressure_mag)*xdiff + (AO.uy(ix,iy,iz)/Pressure_mag)*ydiff + (AO.uz(ix,iy,iz)/Pressure_mag)*zdiff;  % Calculate the unit vector of the US direction

p0 = AO.Pressure(ix,iy,iz);                                                                                                    % Pressure magnitude at this point
phi_a = AO.Phase(ix,iy,iz);                                                                                                    % Phase at this point

Pd_cosj = Pd_cosj + (Tissue.k0*n(n_layer+1)/(Layer.ka*Layer.rho*Layer.va*Layer.va))*dot_a_o*p0*sin(phi_a);                     % Calculate the phase accumulation 
                    
Pd_sinj = Pd_sinj + (Tissue.k0*n(n_layer+1)/(Layer.ka*Layer.rho*Layer.va*Layer.va))*dot_a_o*p0*cos(phi_a);                      % Calculate the phase accumulation 


end