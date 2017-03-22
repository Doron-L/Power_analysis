# Power analysis

This function makes power analysis. It plots the power vs. the effect size for different population sizes. Note that this code doesn't take standard deviation of the measurements into consideration.

Example:  
frac = 87/1328;            % the fraction of the population where you expect to see an effect.  
N = [1e3,5e3,1e4,5e4,1e5]; % Population size vector.  
power_analysis(frac,N)

See power_analysis for more input options.
