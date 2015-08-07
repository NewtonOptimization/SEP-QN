function Ninv_x = quasi_bfgs(S, Q, P, V2, V1, V3, V5, V4, x, rho, hDiag)
% Calculate inv(B+rho)*x
% Using the LBFGS like updating method
% The formula is complex, however it's efficient
%


scale = hDiag + rho;
% scale bigger then run time increase, and convergence slower
% but the step length will remain 1

Ninv_x= x / scale;
t_sx = x' * S;
t_px = x' * P;
t_qx = x' * Q;


for i = 1:size(Q,2)
    s = S(:,i);
    p = P(:, i);
    q = Q(:, i);
    v1 = V1(i);
    v3  = V3(i);
    v5 = V5(i);
    v4 = V4(i);
    
    % v2 precomputed and cached
    v2 = V2(i);
    
    % v6 / v2^2
    v6_v22 = v1 ^ 2 / v2 - 2 * rho * v1 / v2 * v5 + rho * v3 ...
        - rho^2 * v4 + rho^2 / v2 * v5^2;
    
    % this could be computed by MapReduce
    sx = t_sx(i);
    px = t_px(i);
    qx = t_qx(i);
    
    % combination
    s_f = (rho * v5 / v2 * px - v1 / v2 * px + sx - rho * qx ) / v6_v22;
    p_f = ((rho*v5-v1) / v2 * (sx - rho * qx) + (rho^2 * v4 - rho * v3) / v2 * px) / v6_v22;
    q_f = ((rho * v1 - rho^2 * v5) * px / v2 - rho * sx + rho^2 * qx) / v6_v22;
    
    Next_Ninv_x = ( s_f * s +  p_f * p +  q_f * q) ;
    
    Ninv_x = Ninv_x + Next_Ninv_x;    
end


% for i = 1:size(Ninv_S,2)
%     s = S(:,i);
%     v_ny = Ninv_Y(:, i);
%     v_ns = Ninv_S(:, i);
%     ys = YS(i);
%     ss  = SS(i);
%     sny = SNY(i);
%     sns = SNS(i);
%     
%     d = D(i);
%     t_below = -ys ^ 2 / d + 2*rho * sny * ys/d - rho * ss + rho^2 * sns ...
%         - rho^2 * (sny * sny) / d;
%     t_below = - t_below;
%         
%     sx = t_sx(i);
%     v_nyx = t_v_nyx(i);
%     v_nsx = t_v_nsx(i);
%     
%     Fun_Ninv_x_tmp = ( (d^2 * sx  + rho * d * sny * v_nyx - d * ys * v_nyx  - d^2 * rho * v_nsx ) * s ...
%         + ((rho*d*sny-d*ys)* sx + (rho^2 * d * sns - rho * d * ss) * v_nyx + (rho*d*ys - rho^2*d*sny) * v_nsx  )* v_ny ...
%         + ((rho * d * ys - rho^2 * d * sny) * v_nyx - rho*d^2 * sx + rho^2*d^2* v_nsx)* v_ns) / t_below / d^2;
%     
%     Ninv_x = Ninv_x + Fun_Ninv_x_tmp;
%         
% end

end
