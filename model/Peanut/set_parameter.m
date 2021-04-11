
        D        = 3;                   % The dimension of the full system
        d        = 2;                   % The dimension of the slow manifold. Later we can estimate this dimenstion
        chi_p    = 5.991;               % 95% Confidence interval for 2D Gaussian distribution
        K_int    = 100;                 % Number of initial points on the manifold.
        N        = 400000;               % Number of short traj starting from the same initial point
        epsilon  = 0.005;               % the small parameter
        dt       = 5*10^-4;             % Simulation time-step
        T_max    = 0.1;               % Time for short trajctories.
        t0       = 0.1; %0.075;               % The relaxtion time.
        T_one    = 1000;              % Time for single long trajectory to learn ini K charts. 
        threshold = [1, 0.001];                % This is R_max determined by the diameter of the slow manifold when setting the value of R or the maximum possible layor
        connectivity_threshold = 3;   %  two landmarks are connected if they are with in r*sqrt(t0)
        explore_threshold      = 1.05;
        option                 = 1;  % 1 for fast mode prj, 2 for orthogonal prj
        modify                 = 0; % 0 use the intercept(t=0) 1 use the tau_min point as the landmark
        LowerBound = 0.05;       %  Lowerbound of the training window
        UpperBound = 0.10;       %  Upperbound of the training window
        relearn_option     = 0;

        
        
        %parameter for simuation
        Nstep      = 5*10^6;    
        gap        = 10;    
        
        % Parameters of RHS equations
        a1          = 4;
        a2          = 8;
        c1          = 2;
        c2          = 0.5;
        c3          = 0.05;
        c4          = 0.4;
        c5          = 0.05;
        c6          = 0.4; %0.5
        
        % Specify the manifold and gradient
        R               = @(theta, phi)  sqrt(a1+a2*(cos(theta)).^2);
        Rprime          = @(theta, phi) -a2*cos(theta).*sin(theta)./sqrt(a1+a2*cos(theta).^2);
        
        
            % set_up an initial point to start
        theta_int               = rand * pi;
        phi_int                 = rand * 2*pi;
        R_int                   = R(theta_int, phi_int);
        X_int                   = [R_int*sin(theta_int)*cos(phi_int), R_int*sin(theta_int)*sin(phi_int), R_int*cos(theta_int) ];
    
        
        % Specify the drift and diffusion in sph coordinate
        r_RHS           = @(r, theta, phi) -(c1./epsilon)./r.*( r - R(theta,phi) );
        theta_RHS       = @(r, theta, phi)   c3* cos(3*theta)./(r.*sin(theta)) ; % - c3* sin(2*theta + pi/4)./r;
        phi_RHS         = @(r, theta, phi)   c5* sin(phi + theta)./r ; %c5 * sin(phi)./r ;
        sigma_r         = @(r, theta, phi)   c2./(sqrt(epsilon).*r);
        sigma_theta     = @(r, theta, phi)   c4.*sin(theta)./r;
        sigma_phi       = @(r, theta, phi)   c6./r;
        
        b_sph           = @(r, theta, phi) [  r_RHS(r,theta,phi);
                                              theta_RHS(r,theta,phi);
                                              phi_RHS(r,theta,phi) ];
        
        Sigma_sph       = @(r, theta, phi) [  sigma_r(r,theta,phi);
                                             sigma_theta(r,theta,phi);
                                             sigma_phi(r,theta,phi) ];
        

        
        % Pre-load parameters to save some time. 
        drift      = @(X) get_drift( X,epsilon,a1,a2,c1,c2,c3,c4,c5,c6 );
        diffusion  = @(X) get_diffusion( X,epsilon,c2,c4,c6 );

        
  
        
        % Specify the drift and diffusion term on the reduced systems. The drift
        % term has two parts b_stra and b_ito. The latter is the ito term.
        
        b_stra          = @(theta, phi) [( Rprime(theta, phi).*sin(theta)+R(theta, phi).*cos(theta) ).*cos(phi);
                                         ( Rprime(theta, phi).*sin(theta)+R(theta, phi).*cos(theta) ).*sin(phi);
                                         ( Rprime(theta, phi).*cos(theta)-R(theta, phi).*sin(theta) )           ] .* theta_RHS( R(theta,phi), theta, phi ) + ...
                                        [ -R(theta, phi).*sin(theta).*sin(phi);
                                           R(theta, phi).*sin(theta).*cos(phi);
                                           0                                                                    ] .* phi_RHS( R(theta,phi), theta, phi );
        
        H_true          = @(theta, phi) [ (  Rprime(theta, phi).*sin(theta)+R(theta, phi).*cos(theta) ).*cos(phi).*sigma_theta( R(theta,phi), theta, phi ), ...
                                            -R(theta, phi).*sin(theta).*sin(phi).*sigma_phi( R(theta,phi), theta, phi );
                                           ( Rprime(theta, phi).*sin(theta)+R(theta, phi).*cos(theta) ).*sin(phi).*sigma_theta( R(theta,phi), theta, phi ) , ...
                                             R(theta, phi).*sin(theta).*cos(phi).*sigma_phi( R(theta,phi), theta, phi );
                                           ( Rprime(theta, phi).*cos(theta)-R(theta, phi).*sin(theta) ).*sigma_theta( R(theta,phi), theta, phi ), ...
                                             0 ] ;
       
        % this part is calculated by mathematica
        b_ito           = @(theta, phi) [ (1/2).*(a1+a2.*cos(theta).^2).^(-5/2).*cos(phi).*sin(theta).*((-1).*c6.^2.*(    ...
                                           a1+a2.*cos(theta).^2).^2+(-1).*c4.^2.*(a1+a2.*cos(theta).^2).*(a1+4.*a2.*      ...
                                           cos(theta).^2).*sin(theta).^2+a1.*a2.*c4.^2.*sin(theta).^4);
                                          (-1/2).*(a1+a2.*cos(theta).^2).^(-5/2).*sin(theta).*(c6.^2.*(a1+a2.*cos(theta)  ...
                                          .^2).^2+c4.^2.*(a1+a2.*cos(theta).^2).*(a1+4.*a2.*cos(theta).^2).*sin(theta)    ...
                                          .^2+(-1).*a1.*a2.*c4.^2.*sin(theta).^4).*sin(phi);
                                          (-1/4).*c4.^2.*cos(theta).*(a1+a2.*cos(theta).^2).^(-5/2).*(2.*a1.^2+           ...
                                          a2.^2+2.*a2.*(3.*a1+a2).*cos(2.*theta)+a2.^2.*cos(4.*theta)).*sin(theta).^2; ];
        
        
        b_true          = @(X)  btrue(X,b_ito, b_stra);                                                   
        Lambda_true     = @(X)  Lambdatrue(X,H_true);
        R_true          = @(X)  Rtrue(X, R);
        H_true_X        = @(X)  H_trueX(X, H_true);
        
        RHS_parameter   = struct(                ...
            'dt',             dt ,               ...
            'T_max',          T_max ,            ...
            'N',              N,                 ...
            't0',             t0,                ...
            'chi_p',          chi_p,             ...
            'threshold',      threshold,         ...
            'connectivity_threshold',connectivity_threshold, ...
            'D',              D ,                ...
            'd',              d ,                ...
            'modify',         modify,            ...
            'R',              R ,                ...
            'diffusion',      diffusion,         ...
            'drift',          drift,             ...
            'UpperBound',     UpperBound,        ...
            'LowerBound',     LowerBound         ...
            );
        
        
        simulator_one                 = @(Sim_parameter) simEuler_one_traj(RHS_parameter,Sim_parameter);
        simulator_no_par              = @(Sim_parameter) simEuler(RHS_parameter,Sim_parameter);
        simulator_par                 = @(Sim_parameter) simEuler_par(RHS_parameter,Sim_parameter);
        
        Exact_parameter      = struct(           ...
            'H_true_X',       H_true_X,          ...
            'd',              d,                 ...
            'D',              D,                 ...
            'R_true',         R_true,            ...
            'b_true',         b_true,            ...
            'H_true',         H_true,            ...
            'Lambda_true',    Lambda_true        ...
        );
      
        
       chart_sim_parameter      = struct ( ...
                                      'X_int',                   [],   ...
                                      'nearest',                 4,            ...
                                      't0',                      t0,    ...
                                      'dt_s',                    t0,  ...
                                      'D',                       RHS_parameter.D,     ...
                                      'd',                       RHS_parameter.d,     ...
                                      'Nstep',                   Nstep,    ...
                                      'gap',                     gap, ...
                                      'connectivity',            [], ...
                                      'explore_threshold',       explore_threshold,   ...
                                      'N',                       N ...
                                     );   
                                 
        MSM_parameter        = struct(     ...
            'step',               1,       ...
            'N_state',            1*10^6,  ...
            'dt_s',               t0       ...
            );
         
        %parameter for relearning
     
        RHS_parameter_relearn = RHS_parameter; % one can change some RHS if you believe the invariant manifold is not changed.
        relearn_parameter   = struct(                      ...
            'N_relearn',          10*N,                     ...
            'iter',               2,                       ...
            'RHS_parameter',      RHS_parameter_relearn,   ...
            'relative_threshold',  [0.005, 0.02],            ...
            'simulator_par',      @(Sim_parameter) simHeun_par(RHS_parameter_relearn,Sim_parameter) ...
            );
        
       % parameter for MFPT
       N_IC   = 15000;

     
       
% file path       
chart_fileName                     = [datapath,'chart.mat'];
chart_part_fileName                = [datapath,'chart_part.mat'];
TranM_fileName                     = [datapath,'TranM.mat'];
FPT_fileName                       = [datapath,'Peanut_FPT.mat'];





