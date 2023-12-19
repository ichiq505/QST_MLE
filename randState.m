function RHO = randState(nqubit)
    d = 2^nqubit;
    T = 2*vec(rand(d))-1;
    A = diag(T(1:d));
    s = d;
    for i=1:d-1
        COEFs = reshape(T(s:s+2*(d-i)-1),[],2);
        DIAG = complex(COEFs(:,1),COEFs(:,2));
        A = A + diag(DIAG,i);
        s = s + 2*(d-i);
    end
    RHO = A'*A;
    RHO = RHO / real(trace(RHO));
end