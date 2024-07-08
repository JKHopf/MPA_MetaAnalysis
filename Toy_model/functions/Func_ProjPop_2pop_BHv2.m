function [N, R] = Func_ProjPop_2pop_BHv2(scenario, OC, SpParas, Noise, ...
                                      decline, FSq, Ninit, Rtime, Area_vec,...
                                      p, Nsim, pns)
                                  
% Projecting the popualtion over time
% Can be used for both fully fished or spatially managed pops
%
% Jess Hopf
% Aug 2023
%
% Builds on v0 by adding in replicate runs
%
% Inputs:
%   - scenario = 'fished' or 'reserves'
%   - SpParas = single row table of parameters for the species
%   - Ninit = initial population size (matrix size = age-max x nsim)
%   - Rtime = the run time of the projection
%   - Area_vec = area in populations
%   - p = number of populations

% ------------------------------------------------------------------------
% Pre-allocate matrices:
    N = NaN(SpParas.A_max * p, Rtime, Nsim);
    
% ------------------------------------------------------------------------    
% Set up demography

% Age, length, weight relationships:
    % Per capita length (cm) at start of age year (L)
    Lengths = Func_Length((1:SpParas.A_max)',...
                             SpParas.L_inf,SpParas.K,SpParas.A0);

    % Fecundity:
    % Per capita eggs produced           
        % Density indep fecundity (per capita) at length at the start 
        % of age year
        % lengths need to be in mm
        % divide by 2, assuming 50:50 sex ratio
        Fun = Func_fecunds((1:SpParas.A_max)',Lengths*10,...
                             SpParas.c, SpParas.d);

    % Set ages which reproduce
     Fun(1:(SpParas.A_mat-1))= 0; 
      
     
% ------------------------------------------------------------------------
% run over time (inc time-dependent parameters)

% Intial conditions:
N(:,1,:) = Ninit;

% Fishing mortality:
    switch scenario
        case 'fished'
          mf = SpParas.mf;  
        case 'reserve'
            if FSq == "no"
                mf = SpParas.mf;
            elseif FSq == "yes"
                mf = SpParas.mf./(1-Area_vec(1));
            end
        case 'unfished'
          mf = 0;
    end
    
% Survival probabilities:
    % survive natural mortality 
    SrU = exp(-SpParas.m);   
    % survive natural & fishing mort 
    % (inc reallocated effort if reserve scenario)
    SrF = exp(-(SpParas.m+mf));
    
% Proportion caught
%     Df = (mf./(SpParas.m+mf)).*(1-exp(-SpParas.m-mf));

    R = NaN(p,Rtime,Nsim);

for t = 1:Rtime-1 
% Demographic matrices (survival):

    % Ares that are not fished
    RR = diag(repmat(SrU,1,(SpParas.A_max-1)),-1);
    
    % Area that are fished 
    RF = diag([repmat(SrU,1,(SpParas.Ac-1)),...
          repmat(SrF,1,SpParas.A_max-SpParas.Ac)],-1);
    
    % Combine pops
    switch scenario
        case 'fished' 
            B = blkdiag(RF,RF);            
        case 'reserve'
            B = blkdiag(RR,RF);   
        case 'unfished'
            B = blkdiag(RR,RR);   
    end
    
    % Transition matrix
    % repeat for multiple simulaitons
    M = repmat(B,1,1,Nsim); 
%     R = NaN(p,1,Nsim);
    
    % add recruits/larvae
    % how this is added depends on whether its an open or closed pop
    switch OC       
    
        case 'closed'
            % for individual populations
            for i = 1:p    
            % M matrix index values for start/end of population
            Mmin = i*SpParas.A_max-SpParas.A_max+1;
            Mmax = i*SpParas.A_max;
                
            % density (number) of incoming larvae 
            % (same for both per unit area)
            S = sum(repmat(Fun,p,1).*N(:,t,:),1:2); 
        
            % DD p.c. survival - Beverton-Holt function (recruits affect each other)
            R(i,t,:) = SpParas.BHa./(1+(SpParas.BHa.*S)./(SpParas.BHb));
            
            % Multiply recruit survival into projection matrix
            % (duplicated as recruits into one pop come from 2)
            % i.e. p.c. number = fecundity * area * p.c. survival * decline                                  
%             M(Mmin,:,:) = repmat(Fun'.*R(i,:,:),1,p); 
            M(Mmin,:,:) = repmat(Fun'.*Area_vec(i).*R(i,t,:).*decline(t),1,p); 
            end
        
        % multiple matirx by abundance vector
        N(:,t+1,:) = pagemtimes(M,N(:,t,:));      
        
        % ---
    
        case 'open'
            % for individual populations
            for i = 1:p    
            % M matrix index values for start/end of population
            Mmin = i*SpParas.A_max-SpParas.A_max+1;
            Mmax = i*SpParas.A_max;
            
            % DD p.c. survival - Beverton-Holt function
            % (same for both pop under well-mixed assumption) 
            R(i) = SpParas.BHa./(1+(SpParas.BHa.*SpParas.ROpen)./(SpParas.BHb));

            % total number per population (and bulk out zeros)
            Rt(Mmin) = SpParas.ROpen.*Area_vec(i).*R(i);
            Rt(Mmin+1:Mmax) = zeros(SpParas.A_max-1,1);
            end
        
        % multiple matirx by abundance vector
        N(:,t+1,:) = pagemtimes(M,N(:,t,:)); 

        % add recruits 
        N(:,t+1,:) = N(:,t+1,:) + Rt'.*decline(t);

        % ---

    end
    
    % add pink noise
    if Noise{1} == "on"
        switch Noise{2}
            case 'none'
            N(:,t+1,:) = N(:,t+1,:) .* [repmat(pns(t,1,:),SpParas.A_max,1);repmat(pns(t,2,:),SpParas.A_max,1)];
            case 'space'
            N(:,t+1,:) = N(:,t+1,:) .* [repmat(pns(t,1,:),SpParas.A_max,1);repmat(pns(t,1,:),SpParas.A_max,1)];
        end
    end

end





