function [Photon] = LaunchPhoton(R_specular)
%% The function LaunchPhoton  initializes all the variables in Photon structure that photon carries with itself throughout the simulation
 
% INPUTS:

% R_specular: The coefficient of specular reflection 

% OUTPUTS:

% Photon: Structure containing all the variables that photon carries
%%
Photon.x = single(20);
Photon.y = single(20);                                       % Variables to store current Coordinates
Photon.z = single(0);

Photon.x_old = single(0);
Photon.y_old = single(0);                                   % Variables to store previous Coordinates
Photon.z_old = single(0);

Photon.ux_old = single(0);
Photon.uy_old = single(0);                                  % Variables to store previous Direction
Photon.uz_old = single(0);

Photon.ux = single(0);
Photon.uy = single(0);                                      % Variables to store previous Direction
Photon.uz = single(1);

Photon.weight = 1 - R_specular;                     % Reduction of weight due to specular reflection
Photon.dead = false;                                % Initialized Photon is not dead
Photon.s = single(0);                                       % Step size
Photon.s_left = single(0);                                  % Step size i.e. left after the photon travels to the boundary

Photon.Layer = single(1);                                   % Initialize the layer of photon 
                                                    % NOTE: The photon is initialized to layer 1 because I have considered layer 0 as air         
                                                    % and layer (n+1) as air (where n is total number of layers)  
Photon.id = single(0);
                                                    
end
