
% store chart_angle
K = length(chart);
chart_angle       = zeros(d, K);
for i = 1:K
             X_int                 = chart{i}.X_int';
             theta_X_int           = atan2( sqrt( X_int(1)^2 + X_int(2)^2 ), X_int(3) );
             phi_X_int             = mod(atan2(X_int(2), X_int(1)),2*pi);
             chart_angle(:,i)      = [phi_X_int;  theta_X_int ];
end



[V,~]       = eigs(TranM',10);
[~,VD]      = eig(TranM');

if V(1,1) < 0 
     V(:,1)=-V(:,1);
end

angle                        = [pi/6, 5/6*pi];
nearest_index                =  dsearchn(chart_angle',angle);    

if V(nearest_index,2)<0
    V(:,2) = -V(:,2);
end
z2   = V(:,2)';

% periodic extension
for i = 1:K
    
    if chart_angle(1,i)>2*pi-0.2
            temp_phi = chart_angle(1,i)-2*pi;
            temp = [temp_phi;chart_angle(2,i) ];
            chart_angle = [chart_angle, temp];
            z2 = [z2, z2(i)];
    end

    if chart_angle(1,i)<0.2 
         temp_phi = chart_angle(1,i)+2*pi;
         temp = [temp_phi;chart_angle(2,i) ];
         chart_angle = [chart_angle, temp];
         z2 = [z2, z2(i)];

    end

end


x         = chart_angle(1,:);
y         = chart_angle(2,:);
xi        = 0:0.005:2*pi;
yi        = min(y):0.002:max(y);
[xq,yq]   = meshgrid(xi,yi);
zi        = griddata(x,y,z2,xq,yq,'natural');
[c,h]     = contourf(xi,yi,zi,[-0.06:0.02:0.06],'ShowText','on');
s         = contourdata(c);

for i = 1:length(s)
    if s(i).level==0
        phi_zero  = s(i).xdata;
        theta_zero = s(i).ydata;
        break
    end
end


well1 = find(V(:,2)>  0.050); %define cyan state 
well2 = find(V(:,2)< -0.050); %define red  state 

 index1                = find(ismember(nearest_store, well1));
 L_IC                  = length(index1);
 index1                = index1(round(linspace(1, L_IC, N_IC)));
 index2                = find(ismember(nearest_store, well2));
 L_IC                  = length(index2);
 index2                = index2(round(linspace(1, L_IC, N_IC)));
 X_int_store1          = X(index1,:);
 nearest_store1        = nearest_store(index1);
 X_int_store2          = X(index2,:);
 nearest_store2        = nearest_store(index2);
 

well_threshold1 =find(V(:,2)>  0.020);
well_threshold2 =find(V(:,2)< -0.020);


chart_angle       = zeros(d, K);
for i = 1:K
             X_int                 = chart{i}.X_int';
             theta_X_int           = atan2( sqrt( X_int(1)^2 + X_int(2)^2 ), X_int(3) );
             phi_X_int             = mod(atan2(X_int(2), X_int(1)),2*pi);
             chart_angle(:,i)      = [phi_X_int;  theta_X_int ];
end


