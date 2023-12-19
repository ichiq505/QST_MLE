function RHO = noisyState(RHO,error)
    d = length(RHO);
    RHO = (1-error)*RHO + error*(eye(d)/d);
end