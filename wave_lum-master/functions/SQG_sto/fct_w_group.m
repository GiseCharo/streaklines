function w_group_0 = fct_w_group (model, k_wave)
% Compute the wave group veclocity without current
%

w_group_0 = (sqrt(model.physical_constant.g)/2) * ...
    bsxfun( @times, sum( k_wave.^2 , 2) .^(-3/4) , k_wave );
