function [Photon,ppl_dem] = stepInTissue(Tissue,Photon,ppl_dem)
%% The function stepInTissue calculates the step size of photon

% INPUTS:

% Photon: Structure to photon variables.
% Tissue: Structure to tissue variables.
% ppl_dem : Variable to store partial path length 

% OUTPUTS:

% Photon: The updated step-size of photons
% ppl_dem : Variable to store partial path length 

%%
[Layer] = Layer_P(Photon,Tissue);                                          % Call the Layer function to get the current layer of photon
u_t = Layer.mu_a + Layer.mu_s;

if Photon.s_left == 0                                                      % Checks if step size is zero
    rn = rand();                                                           % generate random number
    Photon.s = -log(rn)/(u_t);                                             % Calculate the step size
    ppl_dem(1,Photon.Layer) = ppl_dem(1,Photon.Layer) + Photon.s;           % Store the partial path length with respect to layers
else                                                                       % Checks if step size is not zero
    Photon.s = Photon.s_left/(u_t);                                        % calculate the step size from the left step size
    ppl_dem(1,Photon.Layer) = ppl_dem(1,Photon.Layer) + Photon.s;                                          % Make s_left to zero
    Photon.s_left = single(0);    
end

end
