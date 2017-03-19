function power_analysis(frac,N,varargin)
% (frac,N,effect_size_vec,Effect_size_in_population,alpha,Nsim)

% This function makes power analysis. It plots the power vs. the effect
% size for different population sizes. Note that this code doesn't take standard deviation of the
% measurements into consideration.

% Input:
%------
% frac: the fraction of the population where you expect to see an effect. This can be the fraction that got some treatment
% while the rest is the control, or the fraction of cell with specific mutation, etc.
% N: is the size of the population. A vector with values can be inserted.
% Effect_size_in_population: if one wants to see a change in some existing effect, he should insert here the existing
% effect size in percentage. For example if we want to see if under some condition there is a change in cell proliferation rate,
% which in this case is 1%, insert 0.01. If not, leav3 it in the default value (1).
% alpha: the significance level
% Nsim: number of random binomial simulations,

% Example:
% frac = 87/1328; N = [1e3,5e3,1e4,5e4,1e5]; effect_size_vec = [-100:2:100]/100; Effect_size_in_population = 0.01; alpha = 0.05; Nsim = 1e2;
% power_analysis(frac,N,effect_size_vec,Effect_size_in_population,alpha,Nsim)


    
%% Setting default values
numvarargs = length(varargin);

% set defaults for optional inputs
optargs = {[-100:2:100]/100 1 0.05 1e2};

% now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
effect_size_vec = optargs{1};
[~,Effect_size_in_population, alpha, Nsim] = optargs{:};
%%

if Nsim >= 1e4, fprintf('\n\n Notice, the Nsim value is high, so computation time will be high! \n\n'); end

% Effect size
effect_size_vec = effect_size_vec*Effect_size_in_population; % [%] 

power = zeros(length(N),length(effect_size_vec));
for k = 1:length(N)
    % Make Nsim simulations. In each simulation, we randomly get the number
    % Nk_plus out of N(k) cells using a binomially probability distribution 
    % where the chance for success is frac. 
    Nk_plus = binornd(N(k),frac,[1,Nsim]);
    
    % The number of Nk_minus is simply the complement.
    Nk_minus = N(k) - Nk_plus;
    
    % Simulating the Nk_plus and Nk_minus proliferating cells using the effect
    % size (and not the data).
    Nk_minus_Effect_size = binornd(Nk_minus,Effect_size_in_population);
    
    for j = 1:length(effect_size_vec)
        P_proliferation_plus = Effect_size_in_population + effect_size_vec(j);
        Nk_plus_Effect_size = binornd(Nk_plus,P_proliferation_plus);
        
        Pval = zeros(1,Nsim);
        for i = 1:Nsim
            % The Nk_plus_not_prol = Nk_plus(i)-Nk_plus_prol(i)
            sim_data = [Nk_plus(i)-Nk_plus_Effect_size(i) Nk_plus_Effect_size(i); Nk_minus(i)-Nk_minus_Effect_size(i) Nk_minus_Effect_size(i)];
            % Checking if there is a non random association between Nk_plus and
            % proliferation.
            [~,p] = fishertest(sim_data);
            
            Pval(i) = p;
        end
        % Power is true positive rate, the number of simulations gave
        % P-value more extreme or equal to the significance level (alpha).
        power(k,j) = nnz(Pval<=alpha)/Nsim;
        
    end
end

%% Plotting
figure
plot(smooth(power(1,:),10),'b','LineWidth',2)
hold on
if length(N) > 1
   for i = 2:length(N)
       plot(smooth(power(i,:),10),'LineWidth',2)
   end
end
xlabel('Effect Size','FontSize',15)
ylabel('Power','FontSize',15)
set(gca,'XTick',1:ceil(length(effect_size_vec)/11):length((effect_size_vec)));
set(gca,'XTickLabel',-100:20:100);
set(gca,'XMinorTick','on','YMinorTick','on')
% Make legend   
s1 = [];
for i = 1:length(N)
    s1 = strcat(s1,',''%d''');
end
s1(1) = [];
s2 = sprintf(s1,N);
s3 = sprintf('legend(%s)',s2); eval(s3)


