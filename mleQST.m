function RHO = mleQST(MEAS_DATA)

%     M = repmat(sum(data),4,1);
%     f = data./M;
%     clear('raw_data','M');

    d = uint16(size(MEAS_DATA,2));

%% initialization and setting
    e = 1e-5;

    RHO_pre = eye(d)/d;
    RHO_now = RHO_pre;

    while(1)
        R = evalR(f,RHO_pre);
        E_now = R*RHO_pre*R';
        
        if(cond(RHO_pre,RHO_now) < e)
           break; 
        end

        RHO_pre = RHO_now;
    end
end
function R = evalR(f,RHO_pre)
    % constant variables
    d = length(RHO_pre);
    A = {[1 0; 0 0], [0 0; 0 1], [1 1; 1 1]/2, [1 -1; -1 1]/2, [1 -1i; 1i 1]/2, [1 1i; -1i 1]/2};
    AxAxAxA = cell(d,1);
    for i=1:6
        for j=1:6
            for m=1:6
                for n=1:6
                    AxAxAxA{216*(i-1)+36*(j-1)+6*(m-1)+n} = kron(kron(A{i},A{j}),kron(A{m},A{n}));
                end
            end
        end
    end

    p = zeros(d,1);
    for i=1:d
        p(i) = trace(RHO_pre*AxAxAxA{i}); 
    end

    
    % inner summation
    R0 = zeros(d);
    for m=1:d
        for j=1:36
            for k=1:36
                R0 = R0 + (f(m,j)/p(m,j))*(f(m,k)/p(m,k))*AxAxAxA{j}*E{m}*AxAxAxA{k};
            end
        end
    end
    R0 = inv(sqrtm(R0));
        
    % evaluate R
    R = cell(1,4);
    for n=1:4
        R1 = zeros(4);
        for i=1:36
            R1 = R1 + (f(n,i)/p(n,i)) * R0 * AxAxAxA{i};
        end
        R{n} = R1;
    end
end
function diff = cond(RHO_pre,RHO_now)
    diff = norm(RHO_pre - RHO_now, 'fro');
end