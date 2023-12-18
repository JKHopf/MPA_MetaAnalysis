% MA_Toy_model_v0.m
%
% Jess Hopf
% jess.k.hopf@gmail.com
% July 2023
%
% Purpose:
%   - develop a toy model that demonstrates counterintuative effects when
%     comaparing size, abundance density and biomass density response
%     ratios in MPA analysis
%   - Model is base on California Sheephead (Semicossyphus pulcher)
%
% Model overview:
% 2 population matrix model (MPA/reserve & fished)
% Closed or open population dynamics
% Density-dep model: Beverton-Holt (intra-cohort) recrtuiment
% Fishing with fishery squeeze at time of MPA (CURRENTLY TURNED OFF!)
% Age-structured model

%% Model

% clear workspace
clear 

% add file with relevant functions
addpath('.\functions')

% Variable variables -------------------------------------    

    % number of simulations
    Nsim =  2000; % 200; %  

    % run-times 
    tInit = 200; 
    tPost = 30; 

    % is the population open or closed? 
    OC =  'open'; %   'closed'; %  

    % noise on or off
    % (this is mainly for testing and debugging model)
    Noise{1} = "on"; %  "off"; %  

    % correlated noise?
    % none = pink noise is applied indepdenently to each pop
    % space = pink noise is the same in each pop
    Noise{2} =  'none'; % 'space'; %    

    % declining larval input
    % 1,1 is not declining
    Dec = linspace(1, 1, tInit+tPost);
    % adding a 50% decline over 50 years.
%     Dec(end-49:end) = linspace(1, 0.5, 50);

    % fishery squeeze? 
    FSq =  "no";%   "yes";%


% Parameters ---------------------------------

    % number of populations/patches:
    p = 2; 
    
    % Proportion of area in reserves
    Area = 0.2; %  

    % load species demographic info 
    load SppParameterVals.mat

    % select species
    SpParas = SppParas(SppParas.sp == "Blue rockfish",:);

    % fecundity at length paras (mm):
    % from Dick et al. (2017, Fish Res), Table 6
    SpParas.c = exp(-15.561);
    SpParas.d = 4.816;

    % change fishing pressure 
    % (to match Nickols et al 2019 @ White Rock - medium fishing pressure)
    SpParas.mf = 0.1;

% Noise -------------------------------
    
    % pinknoise
    % this is automatically bounded between [-1,1] though most values fall
    % between [-0.1. 0.1]
    pn = pinknoise(tInit+tPost, Nsim*p);
    
    % scale and then center on 1
    % only scale for non-correlated. If correlated then population
    % fluctuates wildly and often goes extinct, which doesn't match
    % non-correlated as well
    switch Noise{2}
        case 'none'
        pn = pn.*5 + 1;
        case 'space'
        pn = pn.*2.5 + 1;
    end
    
    % reshape
    pns = reshape(pn, size(pn, 1), p, Nsim);


% Model set-up -----------------------------    

    % Model Function:
    DD_func =  str2func('Func_ProjPop_2pop_BHv1'); 
    
    % Closed pop:
        % Beverton-Holt Para values
        % a = slope at/near zero/origin
        % calculuate life-time egg production
        Fun = Func_fecunds((1:SpParas.A_max)',...
              Func_Length((1:SpParas.A_max)',SpParas.L_inf,SpParas.K,SpParas.A0)*10,...
                                 SpParas.c, SpParas.d);
        Fun(1:(SpParas.A_mat-1)) = 0;                      
        LEP = sum(cumprod([1;repmat(exp(-SpParas.m), SpParas.A_max-1, 1)]).*Fun);
        % set so that if population drops below 25% of the unfished LEP then
        % population declines
        SpParas.BHa = 1/(0.25*LEP); % 
        
        % b = max density of settlers
        SpParas.BHb = 1000;    
    
    % Open population: density of recruits coming in
    SpParas.ROpen = 3.2e+11; % 0.2837;
    
    
% Pre-reserv (intial conditions)-----       

    % Fished initial conditions:
    [NPre, NbioPre, R] = DD_func('fished', OC, SpParas, Noise, ...
                        Dec(1:tInit), FSq, ...
                        repmat(10,SpParas.A_max*p,1,Nsim),...
                        tInit, [Area, (1-Area)], p, Nsim, ...
                        pns(1:tInit,:,:));
                
%     figure(100)
%     hold on
%     plot(1:tInit, squeeze(sum(NPre,1)))


% Post-reserve ----
    [NPost, NbioPost, ~] = DD_func('reserve', OC, SpParas, Noise, ...
                        Dec((tInit+1):end), FSq, ...
                        NPre(:,end,:,:),...
                        tPost, [Area, (1-Area)], p, Nsim, ...
                        pns((tInit+1):end,:,:));

%     figure(500)
%     hold on
% %     plot(1:tPost, squeeze(sum(NPost,1)),'k')
%     plot(1:tPost, squeeze(sum(NPost(1:53,:,:),1))./Area)
%     plot(1:tPost, squeeze(sum(NPost(54:106,:,:),1))./(1-Area),'--')
%     legend('R1','R2','F1','F2')
%     
%     figure(302)
%     hold on
%     plot(1:(tInit+tPost), [squeeze(sum(NPre,1));squeeze(sum(NPost,1))])


        
    
    
% Separate populations    
    % index values for start/end of population
    Rind(1) = 1*SpParas.A_max-SpParas.A_max+1;
    Rind(2) = 1*SpParas.A_max;  
    Find(1) = 2*SpParas.A_max-SpParas.A_max+1;
    Find(2) = 2*SpParas.A_max;  
    
    % Totals 
    % Abundance
    RPost = squeeze(sum(NPost(Rind(1):Rind(2),:,:)));
    FPost = squeeze(sum(NPost(Find(1):Find(2),:,:)));
    
    % Density (abundance/unit area)
    Rden = RPost./Area;
    Fden = FPost./(1-Area);
    
    % Biomass density
    Rbioden = squeeze(sum(NbioPost(Rind(1):Rind(2),:,:)./Area));
    Fbioden = squeeze(sum(NbioPost(Find(1):Find(2),:,:)./(1-Area)));

%     Rbioden = squeeze(sum(NbioPre(Rind(1):Rind(2),(tInit-tPost+1):end,:)./Area));
%     Fbioden = squeeze(sum(NbioPre(Find(1):Find(2),(tInit-tPost+1):end,:)./(1-Area)));

    
% 
%     figure(3438)
%     hold on
%     plot(1:tPost, Rden, 'g')
%     plot(1:tPost, Fden, 'b')

%     figure
%     hold on
%     for i = 1:10
%         subplot(10,1,i)
%         plot(1:tPost, Rden(:,i), 'g')
%         hold on
%         plot(1:tPost, Fden(:,i), 'b')
%         hold on
%     end
% 
%     figure
%     hold on
%     plot(1:tPost, mean(Rden./Fden,2))
    

% figure(343)
% hold on
% plot(1:tPost, median(Rden,2), 'g')
% % plot(1:tPost, max(Rden,[],2), '--g')
% % plot(1:tPost, min(Rden,[],2), '--g')
% plot(1:tPost, median(Fden,2), 'b')
% % plot(1:tPost, max(Fden,[],2), '--b')
% % plot(1:tPost, min(Fden,[],2), '--b')

% % before-after
% figure
% hold on
% plot(1:tInit, mean(squeeze(sum(NPre(Rind(1):Rind(2),:,:)./Area)),2), 'g')
% plot(1:tInit, mean(squeeze(sum(NPre(Find(1):Find(2),:,:)./(1-Area))),2), 'b')
% plot((1:tPost)+tInit-1, mean(squeeze(sum(NPost(Rind(1):Rind(2),:,:)./Area)),2), '--g')
% plot((1:tPost)+tInit-1, mean(squeeze(sum(NPost(Find(1):Find(2),:,:)./(1-Area))),2), '--b')
% 
% 
% % general population trend figs
% figure
% hold on
% plot(1:tPost, mean(Rden,2)./mean(Rden(1,:),2), 'g', LineWidth=2)
% plot(1:tPost, mean(Fden,2)./mean(Fden(1,:),2), 'b', LineWidth=2)
% yline(1,'--k')
% ylabel("Density")
% xlabel("Time (yrs)")
% 
% sample example
% figure
% hold on
% plot(1:tPost, Rden, 'g')
% plot(1:tPost, Fden, 'b')


    
    % size
    % length for each age class
    Lengths = Func_Length((1:SpParas.A_max)',...
                           SpParas.L_inf,SpParas.K,SpParas.A0);
    % weight lengths by proporion of pop in each age
    Rsize = squeeze(sum((NPost(Rind(1):Rind(2),:,:)./Area)./sum((NPost(Rind(1):Rind(2),:,:)./Area),1).*Lengths));
    Fsize = squeeze(sum((NPost(Find(1):Find(2),:,:)./(1-Area))./sum((NPost(Find(1):Find(2),:,:)./(1-Area)),1).*Lengths));
        
%         figure(34)
%         subplot(2,1,1)
%         hold on
%         plot(1:tPost, mean(Rsize,1), 'g')
%         subplot(2,1,2)
%         hold on
%         plot(1:tPost, mean(Fsize,1), 'b')
    

%% Response ratios

% Sample the popualtion (i.e. include measurement error )
% affects all years post reserves 
% see function script for more info
%     % aggregation para (degree of clumping, small = more clumping)
% k = 1; 
% NRsamp = Func_MeasureError(NR, meta, k);
% NFsamp = Func_MeasureError(NF, meta, k);
    
% storage order
% I/O, A/B, (AI/BI)/(AO/BO)

% biomass density
RRBD(:,1,:) = reshape(Rbioden./Fbioden, tPost,[],Nsim);
RRBD(:,2,:) = reshape(Rbioden./Rbioden(1,:), tPost,[],Nsim);
RRBD(:,3,:) = reshape(Rbioden./Rbioden(1,:)./(Fbioden./Fbioden(1,:)), tPost,[],Nsim);


% density
RRD(:,1,:) = reshape(Rden./Fden, tPost,[],Nsim);
RRD(:,2,:) = reshape(Rden./Rden(1,:), tPost,[],Nsim);
RRD(:,3,:) = reshape((Rden./Rden(1,:))./(Fden./Fden(1,:)), tPost,[],Nsim);

% size
RRS(:,1,:) = reshape(Rsize./Fsize, tPost,[],Nsim);
RRS(:,2,:) = reshape(Rsize./Rsize(1,:), tPost,[],Nsim);
RRS(:,3,:) = reshape((Rsize./Rsize(1,:))./(Fsize./Fsize(1,:)), tPost,[],Nsim);

    
%     
% % plot
% figure(134)
% hold on
% plot(mean(RRD,3))
% legend('IO','AB','IOAB')
% % plot(RRBD, '--')
% % plot(RRS, '.-')
% ylabel('Response ratio')    
% xlabel('Years')
%     

% Save as table
RR = [reshape(RRBD,[],1);reshape(RRD,[],1);reshape(RRS,[],1)];
time = repmat(1:tPost,1,3*Nsim*3)';
measure = repmat(repelem(["IO","AB","IOAB"],1,tPost),1,Nsim*3)';
rep = repmat(repelem(1:Nsim,tPost*3),1,3)';
var = repelem(["biomass", "density", "size"],1,Nsim*tPost*3)';
Outputs = table(RR, time, measure, rep, var);

% Save csv
% writetable(Outputs, datestr(now, 'yyyy-mm-dd') + "_" + OC + "_FisherySqueeze_" + FSq + "_NoiseCorr_" + Noise{2} + "_decline_" + Dec(end) + "_V0results_5pn2000_BR.csv")




    
%%
figure
hold on
plot(squeeze(RRD(:,1,:)),'k')
plot(squeeze(RRD(:,2,:)),'r')
plot(squeeze(RRD(:,3,:)),'c')
    
    
    
    
   
    