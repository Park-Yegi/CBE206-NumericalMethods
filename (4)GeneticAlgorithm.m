% HOMEWORK 4
% 20160253 Park Yegi

clear all;

mut_rate = 0.2;     % (1) Mutation rate: 20%
N = 10000;          % (2) Number of iteration = 10,000
prob_cross = 0.5;   % (4) Probability of crossover = 50%
chromo_num = 500;   % (5) Size of population = 500
chromos = zeros(chromo_num, 8);
new_chromos = zeros(chromo_num, 8);  
sum_to_250 = 31375; % 250*251/2 = 31375

% Make range of reproduction probabilty for choosing best 250 by ranking
accumulated = 0;
for j = 1:(chromo_num/2)
    accumulated = accumulated + (251-j)/sum_to_250;
    select_prob(j) = accumulated;
end

% Get initial chromosomes
chromos(:, 1:4) = rand(500, 4) * 300;
chromos(:, 5:8) = rand(500, 4) * 10;

% Calculate fitness for initial chromosomes
for j = 1:chromo_num
    fun_val = 0;
    for k = 1:4
        fun_val = fun_val + (chromos(j, k)-50)^2;
    end
    for k = 5:8
        fun_val = fun_val + (chromos(j, k)-5)^2;
    end
    values(j) = fun_val;
end

% Main Iteration
for i = 1:N
    [sorted_val, sorted_idx] = sort(values);
    best_fit = sorted_val(1);
    best_val = chromos(sorted_idx(1), :);

    % Get best fit 250 chromosomes (No mutation)
    % The order is from low fitness to high fitness
    for j = 1:(chromo_num/2)
        new_chromos(j, :) = chromos(sorted_idx(j), :);
    end
    % Also copy the values 
    % in order to reduce calculation time for same chromosomes
    values(1:250) = sorted_val(1:250);
    
    % Get 250 children with mutation
    for j = (chromo_num/2)+1:chromo_num
        beta1 = rand;
        beta2 = rand;
        
        % Select mom and dad among best fit 250 chromosomes
        % according to the select_prob which is made at the beginning
        select_prob1 = rand;
        select_prob2 = rand;
        for k = 1:(chromo_num/2)
            if (select_prob1 < select_prob(k))
                mom = new_chromos(k, :);
                break;
            end
        end
        
        for k = 1:(chromo_num/2)
            if (select_prob2 < select_prob(k))
                dad = new_chromos(k, :);
                break;
            end
        end
        
        % Make new chromosome
        new_chromos(j, 1:4) = mom(1:4)*(1-beta1) + dad(1:4)*beta1;
        new_chromos(j, 5:8) = mom(5:8)*(1-beta2) + dad(5:8)*beta2;
        
        % Mutation
        for k = 1:8
            is_mute = rand;
            if (is_mute < mut_rate)
                if (k<=4)
                    new_chromos(j, k) = rand*300;
                else
                    new_chromos(j, k) = rand*10;
                end
            end
        end
    end
   
    % Update chromosomes
    chromos = new_chromos;
    
    % Calculate fitness (only for children)
    for j = (chromo_num/2)+1:chromo_num
        fun_val = 0;
        for k = 1:4
            fun_val = fun_val + (chromos(j, k)-50)^2;
        end
        for k = 5:8
            fun_val = fun_val + (chromos(j, k)-5)^2;
        end
        values(j) = fun_val;
    end
end

fprintf('the best fitness function is %f and the best chromosome values are e1=%f, e2=%f, e3=%f, e4=%f, s1=%f, s2=%f, s3=%f, s4=%f\n', best_fit, best_val(1), best_val(2), best_val(3), best_val(4), best_val(5), best_val(6), best_val(7), best_val(8));
