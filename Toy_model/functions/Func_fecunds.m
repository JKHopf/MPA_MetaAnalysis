function [Fecunds] = Func_fecunds(Class,Lengths,c,d)
% Fecundity at length at the start of age year (f)
%
% Jess Hopf
%
% Inputs; Class = vector of ages/stages/lengths etc
%        Lengths = vector of lengths       
%        c = length-fecundity parameters 
%        d = length-fecundity parameters 
%
% Outputs; Fecunds = vector of fecundities
%          Class_Fecundity = 'table' of classes and fecunditites

    for i=1:length(Class)
        Fecunds(i,:) = c*(Lengths(i).^d);
    end        
    
end

