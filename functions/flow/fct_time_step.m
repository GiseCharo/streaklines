function  dt = fct_time_step(velocity)
% Determine the optimal time step
%

switch func2str(velocity)
    case {'couetteplanS','couetteplanNS'}
        h=1;
        U0=1;
        nu=1;
        omega=2*nu*((2*pi/h)^2);

        bound_freq = ((2*pi)/omega)/100;
        dx = h/100;
        bound_CFL = dx/U0;
        dt = min([bound_freq bound_CFL]);
        
    case {'DGyreS','DGyreNS'}
        amplitude=0.1;
        omega=pi/5;

        bound_freq = ((2*pi)/omega)/100;
        L = 1;
        dx = L/100;
        vmax = pi*amplitude;
        bound_CFL = dx/vmax;
        dt = min([bound_freq bound_CFL]);
        
    case {'fct_wake','fct_wake_megaRAM'}
        load([ pwd '/data/wakeflow/file_DNS100_inc3d_2017_07_17_1'],'dt');
        
    otherwise
        error('Unknown function')
end
