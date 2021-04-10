function [data, Cov_store, Mean_store ] = simEulernonlintrans(RHS_parameter,Sim_parameter)

diffusion           = RHS_parameter.diffusion;
drift               = RHS_parameter.drift;
nonlin_trans        = RHS_parameter.nonlin_trans;
nonlin_trans_inv    = RHS_parameter.nonlin_trans_inv;
D                   = RHS_parameter.D;
UpperBound          = RHS_parameter.UpperBound;
LowerBound          = RHS_parameter.LowerBound;
X_int               = Sim_parameter.X_int;
T_max               = Sim_parameter.T_max;
dt                  = Sim_parameter.dt;
N                   = Sim_parameter.N;

t_span              = 0:dt:T_max;
tN                  = length(t_span)-1;

Polar = zeros(N,D,tN+1);
X = zeros(N,D,tN+1);
sqrtdt      = sqrt(dt);

r_store     = randn(D,tN*N);
for j=1:N
    Polar_j       = zeros(D,tN+1);
    X_j           = zeros(D,tN+1);
    X_j(:,1)      = X_int';
    Polar_j(:,1)  = nonlin_trans_inv(X_int');
         for i=1:tN
            Current               = Polar_j(:,i);
            cur_r_store           = r_store(:, (j-1)*tN+i );
            drift_current         = drift(Current) * dt;
            diffusion_current     = diffusion(Current) * cur_r_store * sqrtdt;
            
            Next                  = Current +  drift_current + diffusion_current;
            Polar_j(:, i+1)       = Next;
            X_j(:, i+1)           = nonlin_trans(Next);
            
         end
          Polar(j,:,:) = Polar_j;
          X(j,:,:) = X_j;
        
end

Tr_store        = zeros(1, tN+1);
Cov_store       = cell(1,  tN+1);
Mean_store      = cell(1,  tN+1);
for i = 1:tN+1
    Store            = X(:,:,i);
    Cov_store{i}     = cov(Store);
    Mean_store{i}    = mean(Store);
    Tr_store(i)      = trace(Cov_store{i});
end



data         = struct(                          ...
                    'X_int',        X_int,      ...
                    'tN',           tN,         ...
                    'dt',           dt,         ...
                    'T_max',        T_max,      ...
                    'Tr_store',     Tr_store,   ...
                    'LowerBound',   LowerBound, ...
                    'UpperBound',   UpperBound  ...
                 );
end

