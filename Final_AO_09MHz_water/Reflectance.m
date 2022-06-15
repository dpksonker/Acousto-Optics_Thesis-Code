function [Photon,Tissue] = Reflectance(Photon,Tissue,mod_dem,phi_dem,ppl_dem)
%% The function Reflectance records the angle and coordinate in Reflectance matrix and also records weight of photon escaping the first layer

% INPUTS:

% Photon: Structure to photon variables
% Tissue: Structure to tissue variables
% mod_dem : Variable for storing the magnitude of phase modulation in each
%           layers
% phi_dem : Variable for storing the phase of phase modulation in each
%           layers
% ppl_dem : Variable to store partial path length 

% OUTPUTS:

% Photon: Structure to photon variables
% Tissue: Stores the weight of exiting photon and and also the photon
% exiting angle and coordinate
%%
D1 = Tissue.D1;                                                             % distance of detector D1 from source

r_p1  = sqrt((Photon.x - D1)*(Photon.x - D1) + (Photon.y - D1)*(Photon.y - D1));      % Calculate the distance of a photon from the detector

if  r_p1 <= Tissue.Radius_det                                                      % Check if the photon is inside the detector
    Tissue.count_D1 = Tissue.count_D1 + 1; 
    
    Tissue.mod1(Tissue.count_D1,1) = Photon.id;                                      % Save the photon id
    Tissue.mod1(Tissue.count_D1,2) = mod_dem;                                       % Save the modulation magnitude
    Tissue.phi1(Tissue.count_D1,1) = Photon.id;
    Tissue.phi1(Tissue.count_D1,2) = phi_dem;                                        % Save the phase of the phase modulation
    
    Tissue.Det1(Tissue.count_D1,1) =  Photon.id;
    Tissue.Det1(Tissue.count_D1,2) =  Photon.weight;                            % Save the weight of the detected photon
    Tissue.Det1(Tissue.count_D1,3) = ppl_dem;                                   % Save the partial path length 

end

D2 = Tissue.D2;                                                             % distance of detector D2 from source

r_p2  = sqrt((Photon.x - D2)*(Photon.x - D2) + (Photon.y- D2)*(Photon.y- D2));

if  r_p2 <= Tissue.Radius_det 

    Tissue.count_D2 = Tissue.count_D2 + 1; 
    Tissue.mod2(Tissue.count_D2,1) = Photon.id;
    Tissue.mod2(Tissue.count_D2,2) = mod_dem;
    Tissue.phi2(Tissue.count_D2,1) = Photon.id;
    Tissue.phi2(Tissue.count_D2,2) = phi_dem;

    Tissue.Det2(Tissue.count_D2,1) =  Photon.id;
    Tissue.Det2(Tissue.count_D2,2) =  Photon.weight;
    Tissue.Det2(Tissue.count_D2,3) = ppl_dem;
end
% 
D3 = Tissue.D3;                                                              % distance of detector D3 from source

r_p3  = sqrt((Photon.x - D3)*(Photon.x - D3) + (Photon.y - D3)*(Photon.y - D3));

if  r_p3 <= Tissue.Radius_det 

    Tissue.count_D3 = Tissue.count_D3 + 1; 

    Tissue.mod3(Tissue.count_D3,1) =  Photon.id;
    Tissue.mod3(Tissue.count_D3,2) =  mod_dem;
    Tissue.phi3(Tissue.count_D3,1) =  Photon.id;
    Tissue.phi3(Tissue.count_D3,2) =  phi_dem;

    Tissue.Det3(Tissue.count_D3,1) =  Photon.id;
    Tissue.Det3(Tissue.count_D3,2) =  Photon.weight;
    Tissue.Det3(Tissue.count_D3,3) =  ppl_dem;
end



end


