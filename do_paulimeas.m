function [F,P] = do_paulimeas(X,shots)
%% Description
%   X: quantum state to be measured
%   shots: number of shots for sampling, if shots is less then 1 (may be 0?) then
%   it returns theoretically compute a probabilities
%   mode: 0 (sampling) 1(theoretic)

%% code
    num_qubits = log2(length(X));

    F = zeros(1,6^num_qubits);
    P = zeros(1,6^num_qubits);
    E = cell(1,6^num_qubits);
    E = get_pauliPOVM(num_qubits);

    for povm=1:length(E)
        P(povm) = tr(E{povm}*X);
        F(povm) = samp2freq(tr(E{povm}*X),shots);
    end 
    
    F = F / sum(F,'all');
    P = P / sum(P,'all');
end

function freq = samp2freq(prob,shots)
    samples = rand(1,shots);
    freq = nnz(samples < prob) / shots;
end