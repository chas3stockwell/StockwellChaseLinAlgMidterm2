function m = Create_Poisson_problem_A(N)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
    J = -1*eye(N); %identity matrix  N
    Z = zeros(N); %turn this into a for loop. 
    v = ones(N-1, 1) * -1;
    A = diag(v, 1) + diag(v, -1) + 4*eye(N); %N + N-1 + N-1 = 10 = 3N - 2

    m = zeros(N^2);

    for i = 1:N
        r = (i-1)*N+1:i*(N);
        for j = 1:N
            c = (j-1)*N+1:j*(N);
            if i == j
                m(r, c) = A;
            elseif abs(i - j) == 1
                m(r, c) = J;
            else
                m(r, c) = Z;
            end     
        end
    end