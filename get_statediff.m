function diff = get_statediff(X,Y,type)
%% type:
%       - 'fro': Frobenius norm
%       - 'trace' : trace distance
%       - 'fid' : fidelity := ||sqrt(X)*sqrt(Y)||_tr^2

    if strcmp(type,'fro')
        diff = norm(X-Y,'fro');
    elseif strcmp(type,'trace')
        diff = sum(abs(eig(X-Y))) / 2; % trace distance
%     elseif strcmp(type,'fid')
%         sqX=sqrtm(X);
%         sqXYX=sqX*Y*sqX;
%         diff = tr(sqrtm(sqXYX+sqXYX')/2)^2; % fidelity
    else
        diff = norm(X-Y,'fro');
    end
end