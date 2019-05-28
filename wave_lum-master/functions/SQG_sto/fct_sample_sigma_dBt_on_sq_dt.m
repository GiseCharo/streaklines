function sigma_dBt_on_sq_dt = ...
    fct_sample_sigma_dBt_on_sq_dt(model,sigma,coef_modulation)
%% Simulation of sigma dBt

if ~ model.sigma.sto % Deterministic case
    sigma_dBt_on_sq_dt = 0;
else % Stochastic case
    if ~ ( strcmp(model.sigma.type_spectrum,'EOF') || ...
            strcmp(model.sigma.type_spectrum,'Euler_EOF'))
        % Fourier transform of white noise
        dBt_C_on_sq_dt = fft2( randn( ...
            [ model.grid.MX 1 model.sigma.N_ech]));
        % Multiplication by the Fourier transform of the kernel \tilde \sigma
        fft_sigma_dBt_on_sq_dt = bsxfun(@times,sigma,dBt_C_on_sq_dt);
        clear dBt_C_on_sq_dt
        % Homogeneous velocity field
        sigma_dBt_on_sq_dt = real(ifft2(fft_sigma_dBt_on_sq_dt));
        clear fft_sigma_dBt_on_sq_dt
    else
        sigma_dBt_on_sq_dt = sum( sigma .* ...
            randn( [ 1 1 1 N_ech model.sigma.nb_EOF ]) , 5);
    end
    
    %         % Fourier transform of white noise
    %         dBt_C_on_sq_dt = fft2( randn( [ model.grid.MX 1 N_ech]));
    %         % Multiplication by the Fourier transform of the kernel \tilde \sigma
    %         fft_sigma_dBt_on_dt = bsxfun(@times,sigma_on_sq_dt,dBt_C_on_sq_dt);
    %         clear dBt_C_on_sq_dt
    %         % Homogeneous velocity field
    %         sigma_dBt_dt = real(ifft2(fft_sigma_dBt_on_dt));
    %         clear fft_sigma_dBt_on_dt
    
    % Heterogeneous small-scale velocity
    sigma_dBt_on_sq_dt = bsxfun(@times, sqrt(coef_modulation) , ...
        sigma_dBt_on_sq_dt);
    
    if model.sigma.proj_free_div & any(coef_modulation(:) ~= 1)
        % nrj_before_proj_div = mean(sigma_dBt_on_sq_dt(:).^2)
        sigma_dBt_on_sq_dt = fct_proj_free_div(model,sigma_dBt_on_sq_dt);
        % nrj_after_proj_div = mean(sigma_dBt_on_sq_dt(:).^2)
    end
end