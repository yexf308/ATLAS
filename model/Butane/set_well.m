function [X_int_store] = set_well(X,  N_IC, well_threshold)


phi_X                         =        atan2(X(:,6),X(:,4));
index                         =        phi_X < well_threshold(2) & phi_X > well_threshold(1);
X_int_store                   =        X( index ,   :);
[L_IC,~]                      =        size(X_int_store);
X_int_store                   =        X_int_store(round(linspace(1, L_IC, N_IC)),:);

end

