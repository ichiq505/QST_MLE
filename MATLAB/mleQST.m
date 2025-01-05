function [X_final, diffs_k] = mleQST(num_qubits, frequencies, gamma, t_max, alphas, maxiter)
%% Verification & Initialization
%     X: Quantum state estimator
%     F: experimental frequency of click E given X - vector
%     gamma (), tmax(positive-definite value), alphas(< a1 <= a2 <)
%     maxiter: (optional) the maximum number of iterations (equal to max_k)
    
    % default maxiter
    if nargin==5, maxiter = 10000;    end

    % initial settings
    global E F
    th_diff = 1e-6;
    alphas = sort(alphas);
    dim = 2^num_qubits;
    E = get_pauliPOVM(num_qubits);
    F = frequencies;
    k = 0;
    
    diffs_k = [];
%% diluted ML routine
    while 1
        % loop test
        if (k == 0)
            % Initialization of diluted ML routine
            X = eye(dim) / dim;
            t = max([1,t_max]);
            
            disp("-----<< mleQST Information >>-----")
            disp("* num_qubit: " + num_qubits)
            disp("* isVector(F)? " + isvector(F) + ", size of F: " + length(F))
            disp("* gamma: " + gamma)
            disp("* t and t_max: " + t + ", " + t_max)
            disp("* sorted alphas: " + alphas(1) + ", " + alphas(2))
            disp("* max. iteration: " + maxiter)
            disp("* threshold: " + th_diff)
            disp("----------------------------------")
            
%             break
        else % k > 0
            % Test if stop condition holds
            diff_k = get_statediff(X,xX,'fro');
            diffs_k = [diffs_k diff_k];
            if (diff_k - th_diff < 0) || (k > maxiter)
                disp("function stops:")
                disp("* diff_k: " + diff_k)
                disp("* k: " + k)
                break
            end
            clear diff_k
        end
        
        % directional ingredients
        R = get_gradlogML(X);
        RXR = R*X*R;
        trRXR = tr(RXR);
        Db = (R*X + X*R) / 2 - X;
        Dt = RXR / trRXR - X;
        
        % determine next direction subroutine
        while 1
            q = 1 + 2*t + trRXR * (t^2);
            D = (2/q)*Db + (t*trRXR/q)*Dt;
            if ~IsArmijo(X, D, t, gamma)
                % update t in [a1*t, a2*t]
                t = alphas(1)*t + diff(alphas*t)*rand(1);
                continue
            end
            break % if Armijo condition holds, escape from direction determining loop
        end
        
        % update and continue
        k = k+1;
        xX = X;
        X = X + t*D;
    end
    
    X_final = X;
end

function logML = get_logML(X)
%     X: Quantum state estimator
%     F: experimental frequency of click E given X
    global E F
    Y = 0;
    for i=1:length(F)
        Y = Y + F(i) * log(tr(E{i}*X));
    end
    logML = Y;
end

function gradlogML = get_gradlogML(X)
%     X: Quantum state estimator
%     F: experimental frequency of click E given X
%     R(X): gradlogML
    global E F
    Y = 0;
    for i=1:length(F)
        Y = Y + F(i) * E{i} / (tr(E{i}*X));
    end
    gradlogML = Y;
end

function Armijo = IsArmijo(X, D, t, gamma)
%     X: Quantum state estimator
%     R: Gradient of log-likelihood of X
%     t: tentative step-size for dilution
%     D: Direction of the variant of quantum state estimator
%     global F
    R = get_gradlogML(X);
    ArmijoValue = get_logML(X+t*D) - get_logML(X) + gamma * t * tr(R*D);
    Armijo = (ArmijoValue > 0);
end