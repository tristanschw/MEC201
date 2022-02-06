function [sol, its, hist] = myNewton(f, df, x0, TOL, maxit)
g=f(x0);
its=0;
sol=x0;
hist=[];
    while (abs(g)>=TOL) && (its<=maxit)
        hist=[hist,sol];
        sol = sol - g/(df(sol));
        its=its+1;
        g=f(sol);
    end
    hist=[hist,sol];
end