function [Pn_cosi,Pn_sini,Photon,Tissue,AO] = Accumulate_Mod_n(Photon,Tissue,Bound,Pn_cosi,Pn_sini,AO,n)
 %% The function Accumulate_Mod_n accumulates the phase change due to change in refractive index of the medium
 
% INPUTS:

% Photon: Structure to photon variables.
% Bound:  Structure to medium dimensions.
% Tissue: Structure to tissue variables.
% Pn_cosi: Temporary varible for storing the phase accumulations
% Pn_sini: Temporary varible for storing the phase accumulations
% AO :    Structure to loaded variable matrix of Pressure magnitude, Phase
%         , and direction of US
% n :     Tissue layers refractive index matrix 

% OUTPUTS:

% Tissue: Structure to tissue variables.
% Pn_cosi: Temporary varible for storing the phase accumulations
% Pn_sini: Temporary varible for storing the phase accumulations
% Photon: Structure containing all the variables that photon carries
% AO :    Structure to loaded variable matrix of Pressure magnitude, Phase
%         , and direction of US

% % The applied equation is from the article: [Sakadžić, Sava, and Lihong
% V. Wang. "Correlation transfer equation for ultrasound-modulated multiply scattered light." Physical Review E 74.3 (2006): 036618.]
%%  
[Layer] = Layer_P(Photon,Tissue);                                                   % Retrieve the current layer of the photon                                     
ixd = round(abs(Photon.x/Bound.dx));                                                % Calculate the index of photon in x

if ixd == 0

    ix = single(1);

elseif ixd > (Bound.Nx-1)                                                           % Check if the converted index is in the bounds of the total index

    ix = Bound.Nx-1;                                                                % if not then take the last index   

else
    ix = ixd;                                                                       % if yes then take the converted index                                        

end

iyd = round(abs(Photon.y/Bound.dy));                                                % Calculate the index of photon in y

if iyd == 0

    iy = single(1);

elseif iyd > (Bound.Ny-1)

    iy = Bound.Ny-1;

else
    iy = iyd;

end


izd = round(abs(Photon.z/Bound.dz));                                               % Calculate the index of photon in z

if izd == 0

    iz = single(1);
elseif izd > (Bound.Nz-1)

    iz = Bound.Nz-1;

else
    iz = izd;

end

n_layer = Photon.Layer;                                                            % number of layer

p0 = AO.Pressure(ix,iy,iz);                                                        % Pressure magnitude at this point
phi_a = AO.Phase(ix,iy,iz);                                                        % Phase at this point


Pn_cosi = Pn_cosi + (Tissue.k0*n(n_layer+1)*Tissue.eta/(Layer.rho*Layer.va*Layer.va))*Photon.s*p0*cos(phi_a);               % Calculate the phase accumulation 

Pn_sini = Pn_sini - (Tissue.k0*n(n_layer+1)*Tissue.eta/(Layer.rho*Layer.va*Layer.va))*Photon.s*p0*sin(phi_a);               % Calculate the phase accumulation 


end