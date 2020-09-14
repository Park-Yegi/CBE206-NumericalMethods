% HOMEWORK2
% 20160253 ParkYegi

function y = CBE206_hw2_20160253(m, N)
    % Declare variables b
    b1(N) = 1;
    b2(N) = 1;

    %%%% Section #1: Define the banded matrix A %%%%
    for i = 1:N
        for j = 1:N 
            if ((j<i-m)||(j>i+m))
                A1(i, j) = 0;
            else
                A1(i, j) = randn;
            end
        end
    end
    A2 = A1; % copy A1 to A2
    
    
    tic
    %%%% Section #2: original Gaussian elimination %%%%
    % Make lower triangle matrix
    for i = N:-1:2
        for j = i-1:-1:1 
            ratio = A1(j, i)/A1(i, i);
            for k = i:-1:1
                A1(j, k) = A1(j, k) - ratio*A1(i, k);
            end
            b1(j) = b1(j) - ratio*b1(i);
        end
    end
    
    % Back Substitution
    for i = 1:N
        remainder = 0;
        for j = 1:i-1
            remainder = remainder + A1(i, j)*x1(j);
        end
        x1(i) = (b1(i) - remainder)/A1(i, i);
    end
    wtime = toc;
    
    % Print results
    fprintf('original Gaussian elimination\n');
    fprintf('My original program took %f seconds to run\n', wtime);
    
    
    tic
    %%%% Section #3: Optimize Gaussian Elimination %%%%
    % Make lower triangle matrix
    if (m ~= 0)
        for i = N:-1:2
            if (i > m)
                boundary = i-m;
            else
                boundary = 1;
            end

            for j = i-1:-1:boundary
                ratio = A2(j, i)/A2(i, i);
                for k = i:-1:boundary
                    A2(j, k) = A2(j, k) - ratio*A2(i, k);
                end
                b2(j) = b2(j) - ratio*b2(i);
            end
        end
    end
    
    % Back Substitution(same as section #2)
    for i = 1:N
        remainder = 0;
        for j = 1:i-1
            remainder = remainder + A2(i, j)*x2(j);
        end
        x2(i) = (b2(i) - remainder)/A2(i, i);
    end
    wtime = toc;
    
    % Print results
    fprintf('optimal Gaussian elimination\n');
    fprintf('My optimized program took %f seconds to run\n', wtime);
    fprintf('sum1 = %f, sum2 = %f\n', sum(x1), sum(x2));
end
