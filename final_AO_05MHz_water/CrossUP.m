function [Photon,Tissue] = CrossUP(Tissue,Photon,n,mod_dem,phi_dem,ppl_dem)
%% The function CrossUP calculates the fate of photon if it crosses the -ve z direction layer (the -ve z direction is the up direction)
%%  this function calculates if "the photon is reflected" or "is went out of top most layer" or "it travelled to the next layer".  

% INPUTS:

% Tissue: Structure to tissue variables.
% Photon: Structure to photon variables.
% n:      The row vector of refractive indexes of all the layers
% mod_dem : Variable for storing the magnitude of phase modulation in each
%           layers
% phi_dem : Variable for storing the phase of phase modulation in each
%           layers
% ppl_dem : Variable to store partial path length 

% OUTPUTS:

% Photon: The updated direction of photon or photon is dead
% Tissue: saved weight of photon for diffuse reflectance which went out of the medium
%%
layer = Photon.Layer;                                                      % Call the current layer of photon

n_i = n(layer+1);                                                          % Calculate the refractive index of the current layer
n_t = n(layer);                                                            % Calculate the refractive index of the next layer where the photon is going

[Ref] = Rfresnel(Photon,n_i,n_t);                                          % Call the function Rfresnel to calculate the reflection coefficient

ee = rand();                                                               % Generate the random number

if ee <= Ref || Ref == 1                                                   % Check for reflection at boundary

    Photon.uz = -Photon.uz;                                                % if reflected then just change the direction of photon

else                                                                       % Otherwise

    if Photon.Layer == 1                                                   % Check if the photon is at top-most layer                                 

        [Photon,Tissue] = Reflectance(Photon,Tissue,mod_dem,phi_dem,ppl_dem);                % Then call the function to check for diffuse reflectance

        Photon.dead = true;                                                % Make photon dead

    else                                                                   % Otherwise
        uux = Photon.ux;
        uuy = Photon.uy;                                                   % Store the Direction in new variables
        uuz = Photon.uz;

        a_t = acos(abs(uuz));                                              % Calculate the angle of incidence

        Photon.Layer = Photon.Layer - 1;                                   % Reduce the layer of photon 

        Photon.ux = uux*(n_i/n_t);
        Photon.uy = uuy*(n_i/n_t);                                         % Calculate the new direction of photon in another layer it entered
        Photon.uz = sign(uuz)*cos(a_t);

    end

end