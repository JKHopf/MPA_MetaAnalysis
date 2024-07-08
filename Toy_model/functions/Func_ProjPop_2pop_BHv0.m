function [N, Nbio] = Func_ProjPop_2pop_BHv0(scenario, OC, SpParas,...
                                      Ninit, Rtime, Area_vec, p, Nsim)
                                  
% Projecting the popualtion over time
% Can be used for both fully fished or spatially managed pops
%
% Jess Hopf
% Jul 2023
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
    N = NaN(SpParas.A_max * p, Rtime);
    Nbio = N;
    
% ------------------------------------------------------------------------    
% Set up demography

% Age, length, weight relationships:
    % Per capita length (cm) at start of age year (L)
    Lengths = Func_Length((1:SpParas.A_max)',...
                             SpParas.L_inf,SpParas.K,SpParas.A0);

    % Per capita weight at the length at the start of age year (W)
    Weights = Func_Weights(Lengths,SpParas.y,SpParas.z);

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
N(:,1) = repmat(Ninit,1);

% Fishing mortality:
    switch scenario
        case 'fished'
          mf = SpParas.mf;  
        case 'reserve'
          mf = SpParas.mf;% /(1-Area_vec(1));  
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
    Df = (mf./(SpParas.m+mf)).*(1-exp(-SpParas.m-mf));

   
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
    M = B; 
    R = NaN(1,p);
    
    % add recruits/larvae
    % how this is added depends on whether its an open or closed pop
    switch OC       
    
        case 'closed'
            % for individual populations
            for i = 1:p    
            % M matrix index values for start/end of population
            Mmin = i*SpParas.A_max-SpParas.A_max+1;
            Mmax = i*SpParas.A_max;
                
            % density of incoming larvae 
            % (same for both pop under well-mixed assumption) 
            S = sum(repmat(Fun,p,1).*N(:,t));
        
            % DD p.c. survival - Beverton-Holt function (recruits affect each other)
            R(i) = SpParas.BHa./(1+(SpParas.BHa.*S)./(SpParas.BHb));
            
            % Multiply recruits into projection matrix
            % i.e. p.c. number = fecundity * area * p.c. survival                                    
            M(Mmin,:,:,:) = repmat(Fun'.*Area_vec(i).*R(i),1,p); 
            end
        
        % multiple matirx by abundance vector
        N(:,t+1,:) = pagemtimes(M,N(:,t,:));      
    
    
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

            
    
%             % Multiply recruits into projection matrix
%             M(Mmin,:,:,:) = repmat(SpParas.ROpen.*Area_vec(i).*R(i), ...
%                                    1, SpParas.A_max*p); 
    
            end
        
        % multiple matirx by abundance vector
        N(:,t+1,:,:) = pagemtimes(M,N(:,t,:,:)); 

        % add recruits 
        N(:,t+1,:,:) = N(:,t+1,:,:) + Rt';

    end
end


    Nbio = N.*repmat(Weights,p,Rtime);
    
% % Yeild <- OLD to fix before using
%    switch scenario
%         case 'fished' 
%             Y = repmat([zeros(SpParas.Ac-1,size(N,2)); % fished ages only
%                     repmat(Df,SpParas.A_max-SpParas.Ac+1,size(N,2))],4,1).*N; 
%             Ybio = repmat([zeros(SpParas.Ac-1,size(N,2)); % fished ages only
%                     repmat(Df,SpParas.A_max-SpParas.Ac+1,size(N,2))],4,1).*Nbio; 
%         case 'reserve'
%             Y = [zeros(SpParas.A_max*2,size(N,2)); % pops 1 & 2, fished ages only
%                  repmat([zeros(SpParas.Ac-1,size(N,2)); % pops 3 & 4, fished ages only
%                     repmat(Df,SpParas.A_max-SpParas.Ac+1,size(N,2))],2,1)].*N; 
%             Ybio = [zeros(SpParas.A_max*2,size(N,2)); % pops 1 & 2, fished ages only
%                     repmat([zeros(SpParas.Ac-1,size(N,2)); % pops 3 & 4, fished ages only
%                         repmat(Df,SpParas.A_max-SpParas.Ac+1,size(N,2))],2,1)].*Nbio;
%         case 'unfished'
%             Y = 0.*N;
%             Ybio = 0.*Nbio;
%     end



