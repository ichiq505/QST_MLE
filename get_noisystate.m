function RHO = get_noisystate(X,error)
    d = length(X);
    RHO = (1-error)*X + error*(eye(d)/d);
end