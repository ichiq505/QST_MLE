function RHO = get_randstate(num_qubits)
    d = 2^num_qubits;
    t = rand(d);
    T = 2*t(:)-1;
    A = diag(T(1:d));
    s = d;
    for i=1:d-1
        COEFs = reshape(T(s:s+2*(d-i)-1),[],2);
        DIAG = complex(COEFs(:,1),COEFs(:,2));
        A = A + diag(DIAG,i);
        s = s + 2*(d-i);
    end
    RHO = A'*A;
    RHO = RHO / tr(RHO);
end