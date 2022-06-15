function [Layer] = Layer_P(Photon,Tissue)
%% The function Layer_P  initializes the number of layers with their thickness and optical properties
% INPUTS:

% Photon: Structure to photon variables
% Tissue: Structure to tissue variables

% OUTPUTS:

% Layer: Structure to Layer parameters

% NOTE: These optical properties are taken from
% Fang, Xiang, et al. "Effect of scalp hair follicles on NIRS quantification by monte carlo simulation and visible chinese human dataset." 
% IEEE Photonics Journal 10.5 (2018): 1-10.
%%

  f0 = Tissue.f0;
if Photon.Layer == 1                                            % Check if the photon is 1st Layer                       1st LAYER is BRAIN tissue 
    Layer.sz0 =  single(0);                                               % The start of the second layer
    Layer.sz1 =  Tissue.d1;                                    % The thickness of second layer
    Layer.mu_a = single(0.035);                                                    % Absorption coefficient of second layer
    Layer.mu_s = single(8.5);                                                     % Scattering coefficient of second layer
    Layer.g = single(0.9);                                                         % Anisotropy coefficient (g) of first layer
    Layer.rho = single(1000);
    Layer.va = single(1500);
    Layer.ka = single(2*pi*f0/1500);

else
                                                                        % Check if the photon is 1st Layer
    Layer.sz0 = single(0);                                             % Top layer
    Layer.sz1 = single(0);                                        % The thickness "d1" of first layer
    Layer.mu_a = single(0);                                           % Absorption coefficient of first layer
    Layer.mu_s = single(0);                                             % Scattering coefficient of first layer
    Layer.g = single(0);                                                % Anisotropy coefficient (g) of first layer
    Layer.rho = single(0);  
    Layer.va = single(0);  
    Layer.ka = single(0);  
end


end