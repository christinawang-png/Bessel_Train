function reguXguess = regulaRoot(Xguess)
    % Computing the derivitive of a form of regularization by summing up
    % absolute values of neighbor variations.
    % Does not work very well.
    % Xguess: a 3d matrix, stands for the current guessed volume
    %% Regularization formula: R = root(dx^2) + root(dy^2)
    [rx, cx, zx]=size(Xguess);
    reguXguess = zeros(rx,cx,zx);
    for z = 1:zx
        cur = Xguess(:, :, z);
        shiftr_neg = circshift(cur, 1, 1);
        shiftr_neg(1, :) = zeros(1, cx);
        shiftr_pos = circshift(cur, -1, 1);
        shiftr_pos(rx, :) = zeros(1, cx);
        shiftc_neg = circshift(cur, 1, 2);
        shiftc_neg(:, 1) = zeros(rx, 1);
        shiftc_pos = circshift(cur, -1, 2);
        shiftc_pos(:, cx) = zeros(rx, 1);
        %r_neg = zeros(rx, cx); r_pos = zeros(rx, cx); c_neg = zeros(rx, cx); c_pos = zeros(rx, cx);
        r_neg = regulaRootHelper(cur,shiftr_neg);
        r_pos= regulaRootHelper(cur,shiftr_pos);
        c_neg= regulaRootHelper(cur,shiftc_neg);
        c_pos= regulaRootHelper(cur,shiftc_pos);
        reguXguess(:, :, z) = (r_neg + r_pos + c_neg + c_pos).*0.5;
    end  
end