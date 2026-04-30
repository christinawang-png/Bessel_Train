function reguXguess = regula (Xguess)
    % Computing the derivitive of Total Variation Squared Regularization.
    % Easier and Faster to compute than actual Total Variation, but the
    % effects are not that good.
    % Xguess: a 3d matrix, stands for the current guessed volume
    %% Regularization formula: R = dx^2 + dy^2
    TP = [0 -1 0; -1 4 -1; 0 -1 0];
    [rx, cx, zx]=size(Xguess);
    reguXguess = zeros(rx,cx,zx);
    for z = 1:zx 
        correclation = conv2(Xguess(:, :, z), TP, "same");
        reguXguess(:, :, z) = correclation;
    end 
end