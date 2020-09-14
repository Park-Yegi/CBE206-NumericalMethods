% HOMEWORK1
% 20160253 Park Yegi

% Function to calculate nth fibonacci number
% integer -> integer
function my_fib = CBE206_hw1_20160253(num)
    if (num < 3)
        my_fib = 1;
    else
        a = 1;
        b = 1;
        for i=1:(num-2)
            temp = a;
            a = b;
            b = temp + b;
        end
        my_fib = b;
    end
end
