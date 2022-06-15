function [Bound,Tissue,AO,M0,M1] = AO_05MHz_4Layer_BRAIN(Bound,Tissue,n,n_photons,W,AO,M0,M1)
%% The function AO_05MHz_4Layer_BRAIN is the main Monte Carlo simulation function it calculates the total phase modulations due to both AO mechanisms and stores the phase modulations 
%% happened in each layers and also stores the partial path length of each detected photons
    
% INPUTS:

% Bound     : Structure to tissue dimensions
% Tissue    : Structure to tissue variables
% n         : Row matrix of the refractive indexes of all the layers
% n_photons : Total number of photons to be simulated
% W         : threshold for minimum weight of the photon while in
%             simulation
% AO :     Structure to loaded variable matrix of Pressure magnitude, Phase
%         ,and direction of US
% M0:      Variable to accumulate the zeroth order modulation
% M1:      Variable to accumulate the first order modulation

% OUTPUTS:

% Bound   : Structure to tissue dimensions
% Tissue  : Structure to updated RR_dz(1,r) diffuse reflected matrix
% AO :     Structure to loaded variable matrix of Pressure magnitude, Phase
%         ,and direction of US
% M0:      Variable to accumulate the zeroth order modulation
% M1:      Variable to accumulate the first order modulation

%%    
for np=1:n_photons                                                                  % Main loop to simulate the number of photons
    
    Pn_cosi = single(0);
    Pn_sini = single(0);                                                            % Variables for storing phase accumulations at each photon step and scattering                                                 
    Pd_cosj = single(0);
    Pd_sinj = single(0);

    mod_dem = zeros(1,(length(n)-2));                                               % temporary variable for storing the magnitude of phase modulation
    phi_dem = zeros(1,(length(n)-2));                                               % temporary variable for storing the phase of phase modulation
    ppl_dem = zeros(1,(length(n)-2));                                               % temporary variable for storing the partial path length

    r_sp = specularR(n);                                                            % Specular Reflection
    [Photon] = LaunchPhoton(r_sp);                                                  % Launch Photon
     
    count_photon = single(0);
     Photon.id = single(np);                                                        % Variable to count
    
     while (Photon.dead == false)                                                   % Traverse Photon till it is not dead
        count_photon = count_photon + 1;                                            % Count varible to count the number of times the single photon is inside the medium        

        [Photon,ppl_dem] = stepInTissue(Tissue,Photon,ppl_dem);                     % Define step by a random no

        [Photon,hit] = HitBoundary(Photon,Tissue);                                  % Check for the boundary hitting

        if  hit == single(1)                                                        % if hitting

            [Photon] = Hop(Photon);                                                 % Move the photon to the boundary

            if  (Bound.xmin <= Photon.x && Photon.x <= Bound.xmax) && ...           % Make a check if photon is still inside the medium
                    (Bound.ymin <= Photon.y && Photon.y <= Bound.ymax)


                if Photon.uz < single(0)                                            % Check if the photon is at upper boundary

                    [Photon,Tissue] = CrossUP(Tissue,Photon,n,mod_dem,phi_dem,ppl_dem);          % Calculate the fate of photon if its reflected, or out of the medium or traveled to next medium

                else                                                                % Otherwise the photon is crossing the lower boundary

                    [Photon,Tissue] = CrossDOWN(Tissue,Photon,n);                   % Calculate the fate of photon if its reflected, or out of the medium or traveled to next medium
                end

            else
                Photon.dead = true;                                                 % Make photon dead
            end

        else

            [Photon] = Hop(Photon);                                                 % if not hitting the boundary then move the photon

            if  (Bound.xmin <= Photon.x && Photon.x <= Bound.xmax) && ...           % Make a check if photon is still inside the medium
                    (Bound.ymin <= Photon.y && Photon.y <= Bound.ymax)
                
                [Pn_cosi,Pn_sini,Photon,Tissue,AO] = Accumulate_Mod_n(Photon,Tissue,Bound,Pn_cosi,Pn_sini,AO,n);      % Accumulate the phase modulations due to change in refractive index

                [Photon,Tissue,dwa] = Drop(Photon,Tissue);                          % Drop weight as absorbance at a new location
                [Photon] = Spin(Photon,Tissue);                                     % Scatter photon to a new location

                [Pd_cosj,Pd_sinj,Photon,Tissue,AO]       =  Accumulate_ScatterD(Photon,Tissue,Bound,Pd_cosj,Pd_sinj,AO,n);      % Accumulate the phase modulations due to scaterer displacement
                
                [Photon,Tissue,M0,M1,AO,mod_dem,phi_dem] =  AO_Modulation(Photon,Tissue,Pn_cosi,Pn_sini,Pd_cosj,Pd_sinj,Bound,M0,M1,dwa,AO,mod_dem,phi_dem);    % calculate the magnitude and phase of the phase modulation 

            else
                Photon.dead = true;                                                 % Make photon dead                                      

            end

            if  Photon.weight < W  &&  Photon.dead == false || count_photon > 1e4   % Check for photon weight or if it is dead or it has been in a medium > 1e4 
                                                                                    % then make it to face the Roulette

                Photon.dead = true;                                                 % Roulette
            end
        end
    end
end

end

