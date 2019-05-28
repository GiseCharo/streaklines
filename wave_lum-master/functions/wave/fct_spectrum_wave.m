function [spectrum,name_plot,int_epsilon] = fct_spectrum_wave(grid,ft,color)
% Compute the spectrum of a function and superimposed a slope -5/3
%

% Color by default
if nargin < 3
    color='b';
end

% Square modulus of the Fourier transform
ft=abs(ft).^2;
ft= sum(ft,3);

% switch model.dynamics
%     case 'SQG'
%         slope_ref = -5/3;
%     case '2D'
% slope_ref = -3;
%     otherwise
%         error('Unknown type of dynamics');
% end
slope_ref = nan;

% Get parameters
MX=grid.MX;
PX=MX/2;
dX=grid.dX;
if any(size(ft)~=MX)
    error('wrong size');
end
if any( mod(MX,2)~=0)
    error('the number of grid points by axis need to be even');
end

persistent idxref MXref dXref kidx

if ~exist('MXref','var') ||  isempty (MXref) || any(MXref ~= MX) ...
        || any(dXref ~= dX)
    MXref=MX;
    dXref=dX;
    
    % Remove aliasing
    ft(PX(1)+1,:)=0;
    ft(:,PX(2)+1)=0;
    %     ft(PX(1),:)=0;
    %     ft(:,PX(2))=0;
    
    %% Wave vector
    kx=1/(grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
    ky=1/(grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
    kx=2*pi/grid.dX(1)*kx;
    ky=2*pi/grid.dX(2)*ky;
    [kx,ky]=ndgrid(kx,ky);
    k=sqrt(kx.^2+ky.^2);
    k(PX(1)+1,:)=inf;
    k(:,PX(2)+1)=inf;
    k=k(:);
    
    %% Wave number
    M_kappa=min(grid.MX);
    P_kappa= M_kappa/2;
    d_kappa = 2*pi/sqrt(prod(grid.MX.* grid.dX));
    kidx= d_kappa * ( 0:(P_kappa-1) ) ;
    %     M_kappa=min(MX);
    %     P_kappa= M_kappa/2;
    %     d_kappa = max(1./dX);
    %     kidx=1/(M_kappa)* (0:(P_kappa-1)) ;
    %     kidx=2*pi*d_kappa*kidx;
    
    %% Masks associated with the rings of iso wave number
%     d_kappa = kidx(2) - kidx(1);
    kidx_shift = [ kidx(2:end) (kidx(end)+d_kappa) ];
    if M_kappa > 2048
        row_idx = [];
        col_idx = [];
        for i_kidx_local = 1:P_kappa
            kidx_local = kidx(i_kidx_local);
            idx_local = (kidx_local <= k) & (k < kidx_shift(i_kidx_local) ) ;
            idx_local = find( idx_local ) ;
            row_idx = [row_idx ; idx_local];
            col_idx = [col_idx ; i_kidx_local * ones( length(idx_local) ,1 )];
        end
        idx = sparse( row_idx ,col_idx, ...
            ones( length(row_idx) ,1 ) , ...
            prod(grid.MX), P_kappa) ;
    else
        idx = sparse( bsxfun(@le,kidx, k ) );
        idx = idx & sparse( bsxfun(@lt,k, kidx_shift ) );
%         idx = idx & sparse( bsxfun(@lt,k, [ kidx(2:end) (kidx(end)+d_kappa) ] ) );
% %         idx = idx & sparse( bsxfun(@lt,k, (kidx+d_kappa) ) );
    end
    idxref=idx;
    
    %% debug
%     idx_ref=idx;
%     
% %     idx = sparse(prod(grid.MX),P_kappa);
%     clear idx
%     row_idx = [];
%     col_idx = [];
%     for i_kidx_local = 1:P_kappa
%         kidx_local = kidx(i_kidx_local);
%         idx_local = (kidx_local <= k) & (k < kidx_shift(i_kidx_local) ) ;
%         idx_local = find( idx_local ) ;
%         row_idx = [row_idx ; idx_local];
%         col_idx = [col_idx ; i_kidx_local * ones( length(idx_local) ,1 )];
%     end
%     idx = sparse( row_idx ,col_idx, ...
%                   ones( length(row_idx) ,1 ) , ...
%                   prod(grid.MX), P_kappa) ;
%     
%     res = (idx_ref - idx).^2;
%     res = sum(res(:))
%     
    
end

%% Spectrum
% Integration over the rings of iso wave number
spectrum = idxref' * ft(:);

% if strcmp(model.dynamics,'2D')
%     spectrum = spectrum ./(kidx').^2;
% end

% Division by prod(grid.MX) because of the Parseval theorem for
% discrete Fourier transform
% Division by prod(grid.MX) again in order to the integration
% of the spectrum over the wave number yields the energy of the
% buoyancy averaged (not just integrated) over the space
spectrum = 1/prod(grid.MX)^2 * spectrum;

% Division by the wave number step
d_kappa = kidx(2)-kidx(1);
spectrum = spectrum / d_kappa;

% Time integral of the dissipation per scale epsilon(k)
int_epsilon = - cumsum(spectrum) * d_kappa;

%% Plot
idx_not_inf=~(isinf(log10(spectrum(2:end))) ...
    | spectrum(2:end)<1e-4*max(spectrum(2:end)) | isinf(kidx(2:end)'));
idx_not_inf = [false; idx_not_inf];
line1= slope_ref * log10(kidx(2:end))  ;
offset = -1 + mean(  log10(spectrum(idx_not_inf)')  ...
    - line1(idx_not_inf(2:end)));
line1 = line1 + offset;
ref=10.^line1;
loglog(kidx(2:end),ref,'--k');
hold on;
name_plot = loglog(kidx(2:end) , spectrum(2:end) ,color);
ax=axis;
ax(4)=max([spectrum(2:end); ref']);
min_ax= 10 ^(slope_ref * log10(kidx(2)*512/2) + offset) ;
ax(3) = (0.1/(1e-3))^2 * ...
    6e-2*(kidx(2)/kidx(end))*min([max(spectrum); max(ref)']);
% ax(3)=6e-2*(kidx(2)/kidx(end))*min([max(spectrum); max(ref)']);
ax(3) = min( [ax(3) min(ref) min(spectrum(2:end)) min_ax]);
ax(1:2)=kidx(2)*[1 min(grid.MX)/2];
if ax(4)>0
    axis(ax)
end
