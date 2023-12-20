function MEAS_DATA = pauliMeas(RHO)
    d = length(RHO);
    nqubit = log2(d);
    nmeas = 6^nqubit;
    SHOTS = 10000;
    
    pauli = {[1 1;1 1]/2,[1 -1;-1 1]/2,[1 -1i;1i 1]/2,[1 1i;-1i 1]/2,[1 0;0 0],[0 0;0 1]};
    
    %%....%%
    MEAS_DATA = zeros(nmeas,d);
    PROBS = zeros(1,d);
    
    for m=1:nmeas
        MEAS_INDEX = dec2base(m,6,log2(d));
        EFFECT = pauli(uint16(str2double(MEAS_INDEX(1))));
        
        for k=2:nqubit
            EFFECT = kron(EFFECT, pauli(uint16(str2double(MEAS_INDEX(k)))));
        end
                
        outcome = zeros(d);
        for o=1:nqubit
             d
        end
        
        clear EFFECT;
    end
end

function freq = samp2freq(prob,SHOTS)
    samples = rand(SHOTS,1);
    freq = nnz(samples <= prob) / SHOTS;
end