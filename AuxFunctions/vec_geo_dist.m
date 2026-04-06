function v=vec_geo_dist(X)
%%% this function takes a trajectory X, where each column is a coordinate 
%%% for example, for 2-D trajectory, X=[x y], where x and y are column
%%% vector of length equal to the length of the trajectory
%%% the function computes a vector of distances of each point on the
%%% trajectory from the start of the trajectory
dold=0;
for i=1:size(X,1)-1
dnew=norm(X(i,:)-X(i+1,:));
v(i)=dold+dnew; 
dold=v(i);
end

end