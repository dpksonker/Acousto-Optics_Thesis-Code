function [Ref] = Rfresnel(Photon,n1,n2)
%% The function Rfresnel calculates the Fresnel Reflectance

% INPUTS:

% Photon: Structure to photon variables.
% n1: Refractive index of current layer where the photon is present
% n2: Refractive index of next layer where the photon is going

% OUTPUTS:

% Ref: Fresnel reflectance coefficient

%% 
u_z = Photon.uz;                                                           % store the z direction of photon in u_z variable
n_tiss = n1;                                                               % Store the refractive index of current layer 
n_air_t = n2;                                                              % Store the refractive index of next layer 

alpha_ii = pi/2;

alpha_i = acos(abs(u_z));                                                  % the angel of current photon with z-axis

t_emp = (n_air_t/n_tiss);

if t_emp < -single(1) 
    t_emp = -single(0.999); 
end
if t_emp > single(1) 
    t_emp = single(0.999); 
end

critical_a = asin(t_emp);                                       % Calculate the critical angle



temp = (n_tiss/n_air_t)*sin(alpha_i);                                      % temporary variable to store the value


if critical_a < alpha_i && alpha_i < alpha_ii                % Check if the angle of photon is greater than the critical angel
    
    Ref = single(1);                                                 % Total internal reflection happens and make Ref = 1
    
end


if temp > 1 || temp < -1                                    % Check if temp variable is > 1 or < -1
    
    Ref = single(1);

else

alpha_t = asin(temp);
    
    AA = sin(alpha_i - alpha_t)*sin(alpha_i - alpha_t);
    BB = sin(alpha_i + alpha_t)*sin(alpha_i + alpha_t);                     % Fresnel reflection
    CC = tan(alpha_i - alpha_t)*tan(alpha_i - alpha_t);
    DD = tan(alpha_i + alpha_t)*tan(alpha_i + alpha_t);
    
    Ref = single(0.5*((AA/BB) + (CC/DD)));
    
end

end