function [dims, chance_dims, min_dims, max_dims, combs] = patt_agg(patterns, varargin)
%  [dims, random_dims, min_dims, max_dims] = patt_agg(patterns, num_rand_runs, N_rand)  
%
%  Pattern aggregation method (patt agg) compares all combinations 
%    of basispatterns in the same high-d space to see how similar they are.  
%
% INPUTS:
%   patterns: (1 x num_conditions) cell, where
%       patterns{icond} (num_variables x num_patterns)
%       - num_variables (e.g., num_neurons) must be the same 
%           across all conditions
%   num_chance_runs (optional): (1 x 1) scalar, the number of runs to
%       compute the dimensionality expected by chance
%       - default: num_chance_runs = 100
%   N_chance (optional): (1 x 1) scalar, the number of dimensions from
%       which to sample the dimensions expected by chance
%       - must satisfy: max_dim <= N_chance <= num_variables
%       - default: N_chance = num_variables
%
% OUTPUTS:
%   dims: (1 x (2^M-num_conditions-1)), dimensionalities 
%       from effective rank (see below) for each combination 
%       of basis patterns
%   random_dims: (num_chance_runs x (2^M-num_conditions-1)), 
%       dimensionalities that one would expect from
%       chance (randomly drawn from N_chance-dimensional space)
%   min_dims: (1 x (2^M-1)), dimensionalities of best overlap
%   max_dims: (1 x (2^M-1)), dimensionalities if completely orthogonal
%   combs: (1 x (2^M-num_conditions-1)), cell array with the conditions
%       used for each combination.
%       e.g., combs{1} are the conditions whose patterns were compared 
%           for dims{1}
%
% details: 
%   computes the effective of [U_A, U_B, U_C, ...] with threshold 0.5, 
%     where U_A, U_B, U_C, ... are basis patterns for conditions A,B,C,...
%
% example:
%   if patterns consisted of {U_A, U_B, U_C}
%   then dims(1) corresponds to rank([U_A, U_B], 0.5)
%   dims(2) <-- rank([U_A, U_C], 0.5)
%   dims(3) <-- rank([U_B, U_C], 0.5)
%   dims(4) <-- rank([U_A, U_B, U_C], 0.5)
%
% author: Ben Cowley, Dec. 2016
%    ref: Cowley, Smith, Kohn, Yu. "Stimulus-Driven Population Activity 
%           Patterns in Macaque Primary Visual Cortex." PLoS Comp Bio, 2016.

    % perform checks on input
    V = varargin;
    if (checks_on_input_failed(patterns, V))
        error('PAM returned without output.  See above warnings for error.');
    end
    
    % parameters
    rank_thresh = 0.5;  % threshold for rank method
    num_variables = size(patterns{1}, 1); 
    num_conds = length(patterns);
    num_chance_runs = 100;
    N_chance = num_variables;
    
    if (length(varargin) >= 1)
        num_chance_runs = varargin{1};
    end
    if (length(varargin) == 2)
        N_chance = varargin{2};
    end
    
    dims = [];
    chance_dims = [];
    min_dims = [];
    max_dims = [];


    % get the number of patterns for each stimulus
    cond_dim = [];
    for icond = 1:num_conds
        cond_dim(icond) = size(patterns{icond},2);
    end
    
    comb_index = 1;
    combs_cond = [];  % keeps track of which combination of conditions are compared
    for inum_combs = 2:num_conds  % number of elements in combination
        
        combs = nchoosek(1:num_conds, inum_combs); % all possible combinations for chosen number of elements
        
        for icomb = 1:size(combs,1)
            
            combs_cond{comb_index} = combs(icomb,:); % keeps track of combined conditions
            
            % actual patterns
            combined_patterns = [patterns{combs(icomb,:)}];
            dims = [dims rank(combined_patterns, rank_thresh)];
            
            % chance patterns, drawn from N_chance-dimensional space
            for irun = 1:num_chance_runs
                combined_chance_patterns = [];
                for icond = combs(icomb,:)
                    R_cond = orth(randn(N_chance, cond_dim(icond)));
                    combined_chance_patterns = [combined_chance_patterns R_cond];
                end
                chance_dims(irun, comb_index) = rank(combined_chance_patterns, rank_thresh);
            end
            comb_index = comb_index + 1;
            
            % min dims
            min_dims = [min_dims max(cond_dim(combs(icomb,:)))];
            
            % max dims
            max_dims = [max_dims min(sum(cond_dim(combs(icomb,:))), num_variables)];
        end
    end

    combs = combs_cond;

end


function flag = checks_on_input_failed(patterns, V)
% returns true if any check fails

    flag = true;
    
    % patterns all have same number of variables, are not empty
    if (isempty(patterns) || length(patterns) == 1)
        warning('Cell array patterns should have length greater than 1.');
        return;
    end
    
    for icond = 1:length(patterns)
        if (isempty(patterns{icond}))
            warning('Each element of patterns should be non-empty.');
            return;
        end
    end
    
    num_vars = size(patterns{1});
    for icond = 2:length(patterns)
        if (num_vars ~= size(patterns{icond},1))
            warning('Each patterns{icond} should have the same number of variables.');
            return;
        end
    end
    
    % number of patterns no greater than num_variables
    for icond = 1:length(patterns)
        if (num_vars < size(patterns{icond},2))
            warning('Number of basis patterns for each condition must be less than or equal to number of variables.');
            return;
        end
    end
    
    % patterns are orthonormal
    for icond = 1:length(patterns)
        UTU = patterns{icond}' * patterns{icond};
        I = eye(size(UTU));
        if (norm(UTU - I, 'fro') > 1e-5)
            warning('Basis patterns for each condition must be orthonormal.');
            return;
        end
    end
    
    
    % check num_chance_runs and N_chance
    if (length(V) >= 1)
        num_chance_runs = V{1};
        
        if (num_chance_runs <= 0)
            warning('num_chance_runs should be greater than 0.');
            return;
        end
    end
    
    if (length(V) == 2)
        N_chance = V{2};
        if (N_chance > num_vars)
            warning('N_chance should be less than or equal to the number of variables');
            return;
        end
        
        max_num_patterns = 0;
        for icond = 1:length(patterns)
            if (max_num_patterns < size(patterns{icond},2))
                max_num_patterns = size(patterns{icond},2);
            end
        end
        if (N_chance < max_num_patterns)
            warning(sprintf('N_chance should be greater than or equal to the maximum number of patterns %d.', max_num_patterns));
            return;
        end
    end
    
    if (length(V) >= 3)
        warning('Too many arguments were input.');
        return;
    end
        

    flag = false;
end