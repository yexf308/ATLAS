function [X] = simEuler_one_traj(RHS_parameter,Sim_parameter)

% This is Euler Method. Usually for Ito Integral
% X is (#saved states) x dimension

diffusion           = RHS_parameter.diffusion;
drift               = RHS_parameter.drift;
D                   = RHS_parameter.D;
X_int               = Sim_parameter.X_int;
T_max               = Sim_parameter.T_max;
dt                  = Sim_parameter.dt;
gap                 = Sim_parameter.gap;


tN                  = round(T_max/dt);

batch               = 10^8;
X                   = zeros(D,tN/gap+1);
X(:,1)              = X_int';
index_store         = 1;
sqrtdt              = sqrt(dt);
Current             = X_int';


j_total    = floor(tN/batch);
rem        = tN - j_total * batch ;
disp(['Total has ',  num2str(j_total), ' batches.'])

for j = 1:j_total
    disp(['Batch No.', num2str(j)])
    r_store     = randn(D, batch) * sqrtdt;
    for i=1:batch
        drift_current         = drift(Current) * dt;
        diffusion_current     = diffusion(Current) * r_store (:, i );
        Current               = Current +  drift_current + diffusion_current;
        if mod(i,gap)==0
            X( :,index_store+1)               = Current;
            index_store                       = index_store +1;
        end
    end
end

r_store     = randn(D, rem) * sqrtdt;
for i=1:rem
    drift_current         = drift(Current) * dt;
    diffusion_current     = diffusion(Current) * r_store (:, i );
    Current               = Current +  drift_current + diffusion_current;
    if mod(i,gap)==0
        X( :,index_store+1)               = Current;
        index_store                       = index_store +1;
    end
    
end

X = X';

end

