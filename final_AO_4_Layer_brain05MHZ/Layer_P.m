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
if Photon.Layer == 1                                            % Check if the photon is 1st Layer                                                 % Check if the photon is 1st Layer
    Layer.sz0 = single(0);                                                         % Top layer Scalp
    Layer.sz1 = Tissue.d1;                                                 % The thickness "d1" of first layer
    Layer.mu_a = single(0.0375);                                                    % Absorption coefficient of first layer
    Layer.mu_s = single(9.5);                                                       % Scattering coefficient of first layer
    Layer.g = single(0.8);                                                         % Anisotropy coefficient (g) of first layer
    Layer.rho = single(1116);
    Layer.va = single(1537);
    Layer.ka = single(2*pi*f0/1537);

elseif Photon.Layer == 2                                                    % layer Skull

    Layer.sz0 = Tissue.d1;
    Layer.sz1 =  Tissue.d1 + Tissue.d2;
    Layer.mu_a = single(0.028);
    Layer.mu_s = single(27);
    Layer.g = single(0.9);
    Layer.rho = single(1796);
    Layer.va = single(2652);
    Layer.ka = single(2*pi*f0/2652);
    %
elseif Photon.Layer == 3                                                    % layer CSF

    Layer.sz0 = Tissue.d1 + Tissue.d2;
    Layer.sz1 = Tissue.d1 + Tissue.d2 + Tissue.d3;
    Layer.mu_a = single(0.00433);
    Layer.mu_s = single(0.0158);
    Layer.g = single(0.9);
    Layer.rho = single(1000);
    Layer.va = single(1500);
    Layer.ka = single(2*pi*f0/1500);
else
    Layer.sz0 = Tissue.d1 + Tissue.d2 + Tissue.d3;                                                % The start of the second layer
    Layer.sz1 =  Tissue.d1 + Tissue.d2 + Tissue.d3 + Tissue.d4;                                    % The thickness of second layer
    Layer.mu_a = single(0.035);                                                    % Absorption coefficient of second layer
    Layer.mu_s = single(8.5);                                                     % Scattering coefficient of second layer
    Layer.g = single(0.9);                                                         % Anisotropy coefficient (g) of first layer
    Layer.rho = single(1000);
    Layer.va = single(1500);
    Layer.ka = single(2*pi*f0/1500);

end