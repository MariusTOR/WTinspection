function [convergence] = newton(equation,dequation,initial,tol)
% Newton-Raphson solver for an anonymous function
% Coded by Br. Marius Strom, TOR
    convergence = initial;
    delta = 1;
    while(delta>tol)
        delta = convergence;
        convergence = convergence - (equation(convergence) - convergence) ./ dequation(convergence);
        delta = max(abs(convergence - delta));
    end
end

