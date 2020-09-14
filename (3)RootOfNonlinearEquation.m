% HOMEWORK #3
% 20160253 Park Yegi

clear all;

% Declare variables needed
syms x y;
num_sol = 0; % This variable is for counting the number of the solutions found.
epsilon = 1e-5;
x_sol = NaN;
y_sol = NaN;

% Declare two functions
f1(x, y) = x^2 + x*y -10;
f2(x, y) = y + 3*x*(y^2) - 57;

% Calculate partial derivatives
df1_dx(x, y) = diff(f1, x);
df1_dy(x, y) = diff(f1, y);
df2_dx(x, y) = diff(f2, x);
df2_dy(x, y) = diff(f2, y);


% main loop for finding solutions
while true
    % initimalize variables.
    num_iter = 1; % This variable is for counting the number of iterations for certain guess.
    norm = 1e-2;

    % Get random number between [-5, 5]
    x_prev = 10.*rand -5;
    y_prev = 10.*rand -5;
    
    
    % loop for solving each system of linear equations
    while (norm > epsilon && num_iter <= 100)
        % Create Jacobian matrix and b vector
        Jacobian = zeros(2, 2);
        Jacobian(1, 1) = double(df1_dx(x_prev, y_prev));
        Jacobian(1, 2) = double(df1_dy(x_prev, y_prev));
        Jacobian(2, 1) = double(df2_dx(x_prev, y_prev));
        Jacobian(2, 2) = double(df2_dy(x_prev, y_prev));
        b_vector(1) = double(f1(x_prev, y_prev)) * (-1);
        b_vector(2) = double(f2(x_prev, y_prev)) * (-1);
        
        % Gaussian Elimination
        % Make upper triangular matrix
        ratio = Jacobian(2, 1)/Jacobian(1, 1);
        Jacobian(2, 1) = Jacobian(2, 1) - ratio*Jacobian(1, 1);
        Jacobian(2, 2) = Jacobian(2, 2) - ratio*Jacobian(1, 2);
        b_vector(2) = b_vector(2) - ratio*b_vector(1);
        % Back substitution
        y_next = b_vector(2) / Jacobian(2, 2);
        x_next = (b_vector(1) - y_next*Jacobian(1, 2)) / Jacobian(1, 1);
        y_next = y_next + y_prev;
        x_next = x_next + x_prev;
        
        % Advance
        norm = sqrt(abs(x_prev-x_next).^2 + abs(y_prev-y_next).^2);
        x_prev = x_next;
        y_prev = y_next;
        num_iter = num_iter +1;
    end


    % if finding the solutions, print them
    if (num_iter <= 100)
        num_sol = num_sol +1;
        if (num_sol == 1)
            fprintf('1st solution: %f, %f\n', x_prev, y_prev);
        elseif (num_sol == 2)
            fprintf('2nd solution: %f, %f\n', x_prev, y_prev);
        elseif (num_sol == 3)
            fprintf('3rd solution: %f, %f\n', x_prev, y_prev);
        else
            fprintf('%dth solution: %f, %f\n', num_sol, x_prev, y_prev);
        end
        
        % When we find the second solution.
        if (~isnan(x_sol) && round(x_prev) ~= round(x_sol))
            break;
        end

        x_sol = x_prev;
        y_sol = y_prev;
    end
end
