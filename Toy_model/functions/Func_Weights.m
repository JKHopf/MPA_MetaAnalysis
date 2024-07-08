function [Weights] = Func_Weights(Lengths,y,z)

% Weights (g) at start of age year (W)
%
% Jess Hopf
%
% Inputs; Age = vector of ages
%         Lenghts = vector of lengths 
%         y = length-weight parameter
%         z = length-weight parameter
% Outputs; Weights = vector of weights
%          Age_Weights = 'table' of age and weights


     for i=1:length(Lengths)
        Weights(i,:) = y*(Lengths(i)).^z;
     end

end
