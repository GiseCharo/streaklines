function model = fct_set_forcing(model)
% set the parameter of the forcing
%

% if isfield(model.advection, 'forcing') && model.advection.forcing.bool
%     %     Ly = model.grid.MX(2) * model.grid.dX(2);
%     %     [~,Y]=ndgrid(model.grid.x,model.grid.y);
%     %     Vy = model.advection.forcing.ampli_forcing * model.odg_b * ...
%     %         sin( 2 * pi * model.advection.forcing.freq_f * Y / Ly);
%     %     clear Ly X Y
%
%
%     model.advection.forcing.Lx = model.grid.MX(1) * model.grid.dX(1);
%     model.advection.forcing.Ly = model.grid.MX(2) * model.grid.dX(2);
%     [X,Y]=ndgrid(model.grid.x,model.grid.y);
%
%     switch model.dynamics
%         case 'SQG'
%             U_caract = ...
%                 model.odg_b / model.physical_constant.buoyancy_freq_N;
%         case '2D'
%             U_caract = ...
%                 model.odg_b * model.advection.forcing.Ly;
%             U_caract = U_caract /20;
%             U_caract = U_caract /4;
%         otherwise
%             error('Unknown type of dynamics');
%     end
%     %     U_caract = U_caract /10;
%     U_caract = U_caract /3;
%
%     %     model.advection.U_caract = sqrt(mean(w(:).^2));
%     model.advection.forcing.on_T =  U_caract / model.advection.forcing.Ly;
%
%     %     model.advection.forcing.on_T = model.advection.forcing.on_T /10;
%
%     model.advection.forcing.ampli_forcing = ...
%         model.advection.forcing.ampli_forcing * model.odg_b ;
%     %     model.advection.forcing.ampli_forcing = ...
%     %         model.advection.forcing.ampli_forcing * model.odg_b * ...
%     %         model.advection.forcing.on_T;
%
%     %     ampli_scale = sqrt( ...
%     %         ( sum(model.advection.forcing.freq_f.^2)/2 ) ...
%     %             ^(model.sigma.slope_sigma) ...
%     %             );
%     ampli_scale = 1;
%
%     switch model.advection.forcing.forcing_type
%         case 'Kolmogorov'
%             F = ampli_scale * ...
%                 model.advection.forcing.ampli_forcing * ...
%                 cos( 2 * pi / model.advection.forcing.Lx ...
%                 * model.advection.forcing.freq_f(1) * X ...
%                 + 2 * pi / model.advection.forcing.Ly ...
%                 * model.advection.forcing.freq_f(2) * Y );
%             F = F / model.advection.forcing.on_T;
%         case 'Spring'
%             F = ampli_scale * ...
%                 model.advection.forcing.ampli_forcing ...
%                 * sin( 2 * pi / model.advection.forcing.Lx ...
%                 * model.advection.forcing.freq_f(1) * X ) ...
%                 .* sin( 2 * pi / model.advection.forcing.Ly ...
%                 * model.advection.forcing.freq_f(2) * Y );
%             %     F = model.advection.forcing.ampli_forcing * ...
%             %         sin( 2 * pi * model.advection.forcing.freq_f * Y / ...
%             %         model.advection.forcing.Ly);
%     end
%
%
%
%     model.advection.forcing.F = fft2(F);
%
%     if strcmp(model.type_data, 'Zero')
%         fft_w = SQG_large_UQ(model,  ...
%             model.odg_b / model.advection.forcing.ampli_forcing ...
%             * model.advection.forcing.F);
%         %   w(:,:,1) = U_caract / model.advection.forcing.ampli_forcing * F;
%         %         w(:,:,2) = 0;
%         %         fft_w = fft2(w);
%         if strcmp( model.advection.forcing.forcing_type,'Kolmogorov')
%             fft_w = fft_w * model.advection.forcing.on_T;
%         end
%         w = real(ifft2(fft_w));
%     end
%
%     %
%     %     figure;imagesc(model.grid.x,model.grid.y,model.advection.forcing.F');
%     %     axis xy; axis equal;
%     %
%
%     %     model.advection.forcing.F = fft2(F);
%     clear Lx Ly X Y on_T U_caract ampli_scale
% end
if isfield(model.advection, 'forcing') && model.advection.forcing.bool
    %     Ly = model.grid.MX(2) * model.grid.dX(2);
    %     [~,Y_forcing]=ndgrid(model.grid.x,model.grid.y);
    %     Vy = model.advection.forcing.ampli_forcing * model.odg_b * ...
    %         sin( 2 * pi * model.advection.forcing.freq_f * Y_forcing / Ly);
    %     clear Ly X_forcing Y_forcing
    
    
    model.advection.forcing.Lx = model.grid.MX(1) * model.grid.dX(1);
    model.advection.forcing.Ly = model.grid.MX(2) * model.grid.dX(2);
    [X_forcing,Y_forcing]=ndgrid(model.grid.x,model.grid.y);
    
    
    
    switch model.dynamics
        case 'SQG'
            U_caract = ...
                model.odg_b / model.physical_constant.buoyancy_freq_N;
            U_caract = U_caract /5;
        case '2D'
            U_caract = ...
                model.odg_b * model.advection.forcing.Ly;
            U_caract = U_caract /20;
            U_caract = U_caract /4;
        otherwise
            error('Unknown type of dynamics');
    end
    %     U_caract = U_caract /10;
    U_caract = U_caract /3;
    
    %     model.advection.U_caract = sqrt(mean(w(:).^2));
    model.advection.forcing.on_T =  U_caract / model.advection.forcing.Ly;
    
    %     model.advection.forcing.on_T = model.advection.forcing.on_T /10;
    
    model.advection.forcing.ampli_forcing = ...
        model.advection.forcing.ampli_forcing * model.odg_b ;
    %     model.advection.forcing.ampli_forcing = ...
    %         model.advection.forcing.ampli_forcing * model.odg_b * ...
    %         model.advection.forcing.on_T;
    
    %     ampli_scale = sqrt( ...
    %         ( sum(model.advection.forcing.freq_f.^2)/2 ) ...
    %             ^(model.sigma.slope_sigma) ...
    %             );
    ampli_scale = 1;
    
    switch model.advection.forcing.forcing_type
        case 'Kolmogorov'
            F = ampli_scale * ...
                model.advection.forcing.ampli_forcing * ...
                cos( 2 * pi / model.advection.forcing.Lx ...
                * model.advection.forcing.freq_f(1) * X_forcing ...
                + 2 * pi / model.advection.forcing.Ly ...
                * model.advection.forcing.freq_f(2) * Y_forcing );
            F = F / model.advection.forcing.on_T;
        case {'Spring','Hetero_Spring'}
            F = ampli_scale * ...
                model.advection.forcing.ampli_forcing ...
                * sin( 2 * pi / model.advection.forcing.Lx ...
                * model.advection.forcing.freq_f(1) * X_forcing ) ...
                .* sin( 2 * pi / model.advection.forcing.Ly ...
                * model.advection.forcing.freq_f(2) * Y_forcing );
            %     F = model.advection.forcing.ampli_forcing * ...
            %         sin( 2 * pi * model.advection.forcing.freq_f * Y_forcing / ...
            %         model.advection.forcing.Ly);
    end
    
    
    if strcmp(model.advection.forcing.forcing_type, 'Hetero_Spring')
        model.advection.forcing.F = F;
    else
        model.advection.forcing.F = fft2(F);
    end
    
    if strcmp(model.type_data, 'Zero')
        fft_w = SQG_large_UQ(model,  ...
            model.odg_b / model.advection.forcing.ampli_forcing ...
            * model.advection.forcing.F);
        %         w(:,:,1) = U_caract / model.advection.forcing.ampli_forcing * F;
        %         w(:,:,2) = 0;
        %         fft_w = fft2(w);
        if strcmp( model.advection.forcing.forcing_type,'Kolmogorov')
            fft_w = fft_w * model.advection.forcing.on_T;
        end
        w = real(ifft2(fft_w));
    end
    
    %
    %     figure;imagesc(model.grid.x,model.grid.y,model.advection.forcing.F');
    %     axis xy; axis equal;
    %
    
    %     model.advection.forcing.F = fft2(F);
    clear Lx Ly X_forcing Y_forcing on_T U_caract ampli_scale
end