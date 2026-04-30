function reguXguess = regulaTV(Xguess)
    % Computing the derivitive of Total Variation regularization
    % Xguess: a 3d matrix, stands for the current guessed volume
    %% Regularization formula: R = root(dx^2 + dy^2)
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

        r_neg = zeros(rx, cx); r_pos = zeros(rx, cx); c_neg = zeros(rx, cx); c_pos = zeros(rx, cx);

        r_pos(1:rx-1, :) = cur(1:rx-1, :) - shiftr_pos(1:rx-1, :);
        r_neg(2:rx, :)  = cur(2:rx, :) - shiftr_neg(2:rx, :);
        c_pos(:, 1:cx-1) = cur(:, 1:cx-1) - shiftc_pos(:, 1:cx-1);
        c_neg(:, 2:cx) = cur(:, 2:cx)  - shiftc_neg(:, 2:cx);

        TVmat = sqrt(r_pos.^2 + c_pos.^2);
        TVmat_shiftr_neg = circshift(TVmat, 1, 1);
        TVmat_shiftr_neg(1, :) = zeros(1, cx);
        TVmat_shiftc_neg = circshift(TVmat, 1, 2);
        TVmat_shiftc_neg(:, 1) = zeros(rx, 1);

        result_1 = (r_pos + c_pos)./ TVmat;
        result_2 = c_neg ./ TVmat_shiftc_neg;
        result_3 = r_neg ./ TVmat_shiftr_neg;
        result_1(~isfinite(result_1))=0;
        result_2(~isfinite(result_2))=0;
        result_3(~isfinite(result_3))=0;

        result = result_1 + result_2 + result_3;
        reguXguess(:, :, z) = result;
    end
end