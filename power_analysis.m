function power_analysis(N,varargin)
% power_analysis(N,effect_size,Effect_size_in_population,Treatment_sample_frac,alpha,Nsim,Treatment_effect_frac)

% This function makes power analysis. It plots the power vs. the effect
% size for different population sizes. Note that this code doesn't take standard deviation of the
% measurements into consideration.

% Input:
%------
% N: is the size of the population. A vector with values can be inserted.
%
% effect_size: a vector of possible effect size values.
% 
% Effect_size_in_population: if one wants to see a change in some existing 
% effect, he should insert here the existing effect size in percentage. For 
% example if we want to see if under some condition there is a change in 
% cell proliferation rate, which in this case is 1%, insert 0.01. If not, 
% leave it in the default value (0). 
% 
% Treatment_sample_frac: the treatment sample fraction of the total 
% sample (treatment+control). The treatment sample can be a sample that some treatment was 
% applied on it or a sample that is different by the control whether it is 
% a mutation. The default value is 1/2 (which means the treatment and control 
% group have the same size.
%
% alpha: the significance level
%
% Nsim: number of random binomial simulations
%
% Treatment_effect_frac: is the fraction within the treatment sample on
% which we expect to see an effect size. The default is all the treatment
% sample (1).

% Example:
% Treatment_sample_frac = 87/1328; N = [1e3,5e3,1e4,5e4,1e5]; effect_size = [-100:2:100]/100; Effect_size_in_population = 0.01; alpha = 0.05; Nsim = 1e2;
% power_analysis(N,effect_size,Effect_size_in_population,Treatment_sample_frac,alpha,Nsim)


    
%% Setting default values
numvarargs = length(varargin);

% set defaults for optional inputs
optargs = {[-100:2:100]/100 0 0.5 0.05 1e2 1};

% now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
effect_size = optargs{1};
[~,Effect_size_in_population, Treatment_sample_frac, alpha, Nsim, Treatment_effect_frac] = optargs{:};
%%

if Nsim >= 1e4, fprintf('\n\n Notice, the Nsim value is high, so computation time will be high! \n\n'); end

% Effect size
% If this is an effect size of some feature/effect that already exists in
% all the population (Effect_size_in_population ~= 0), we take it into
% consideration.
if Effect_size_in_population ~= 0
    effect_size = effect_size*Effect_size_in_population; % [%]
end
% If the effect size doesn't apply on all the treatment sample, we take it
% into consideration.
effect_size = effect_size*Treatment_effect_frac;

power = zeros(length(N),length(effect_size));
for k = 1:length(N)
    % Make Nsim simulations. In each simulation, we randomly get the number
    % Nk_plus out of N(k) cells using a binomially probability distribution 
    % where the chance for success is Treatment_sample_frac. 
    Nk_plus = binornd(N(k),Treatment_sample_frac,[1,Nsim]);
    
    % The number of Nk_minus is simply the complement.
    Nk_minus = N(k) - Nk_plus;
    
    % Simulating the Nk_plus and Nk_minus proliferating cells using the effect
    % size (and not the data).
    Nk_minus_Effect_size = binornd(Nk_minus,Effect_size_in_population);
    
    for j = 1:length(effect_size)
        P_effect_size_plus = Effect_size_in_population + effect_size(j); % Taking into consideration the effect size in the treatment 
                                                                             % sample and the effect size in all population. 
        Nk_plus_Effect_size = binornd(Nk_plus,P_effect_size_plus);
        
        Pval = zeros(1,Nsim);
        for i = 1:Nsim
            % The Nk_plus_not_prol = Nk_plus(i)-Nk_plus_prol(i)
            sim_data = [Nk_plus(i)-Nk_plus_Effect_size(i) Nk_plus_Effect_size(i); Nk_minus(i)-Nk_minus_Effect_size(i) Nk_minus_Effect_size(i)];
            % Checking if there is a non random association between Nk_plus
            % and the effect. For example, 
%                          Flu    NoFlu
%                          ___    _____
% 
%                 NoShot    3      6    
%                 Shot      1      7    
% 
%             Use Fisher's exact test to determine if there is a nonrandom association between receiving a flu shot and getting the flu.
%sim_data
            [~,p] = fishertest(sim_data);
            
            Pval(i) = p;
        end
        % Power is true positive rate, the number of simulations gave
        % P-value more extreme or equal to the significance level (alpha).
        power(k,j) = nnz(Pval<=alpha)/Nsim;
        
    end
end

%% Plotting
% Getting back the effect size that apply on the effect and items we expect
% it to apply.
effect_size = effect_size/Effect_size_in_population/Treatment_effect_frac*1e2;
figure
plot(effect_size,smooth(power(1,:),10),'LineWidth',2)
hold on
if length(N) > 1
   for i = 2:length(N)
       plot(effect_size,smooth(power(i,:),10),'LineWidth',2)
   end
end
xlabel('Effect Size [%]','FontSize',15)
ylabel('Power','FontSize',15)
set(gca,'XMinorTick','on','YMinorTick','on')
% Make legend   
s1 = [];
for i = 1:length(N)
    s1 = strcat(s1,',''%d''');
end
s1(1) = [];
s2 = sprintf(s1,N);
s3 = sprintf('legend(%s)',s2); eval(s3)


