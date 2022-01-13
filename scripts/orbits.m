
addpath(genpath("."));

[T,X] = obtain_2D_motion([6791; 0; 0; 7.66*1.3], [0; 5600*3]);


figure
plot(X(:,1), X(:,2))
daspect([1 1 1])


[T,X] = obtain_3D_motion([6791; 0; 0; 0; 7.66*1.3*cos(0.3); 7.66*1.3*sin(0.3)], [0; 5600*3]);

%Hi

figure
plot3(X(:,1), X(:,2), X(:,3))
daspect([1 1 1])