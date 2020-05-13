using Plots
#make the spring node coordinates
D = 20e-3;
windings = 1;
p = 5e-3;
n = 10;
t = LinRange(0,windings,n);
x = p*t;
y = 0.5*D*cos.(2*pi*t);
z = 0.5*D*sin.(2*pi*t);
plot(x,y,z)

#define topology
xyz = [x;y;z];
numel[:,1] = 1:1:n-1;
numel[:,2] = 2:1:n;

#start main loop to calculate stiffness matrix
for i=1:n-1
    i1 = numel[i,1];
    i2 = numel[i,2];
    #calclulate directions
    x1 = x[i1]; x2 = x[i2]; dx = x2-x1;
    y1 = y[i1]; y2 = y[i2]; dy = y2-y1;
    z1 = z[i1]; z2 = z[i2]; dz = z2-z1;
    L = sqrt(dx^2+dy^2+dz^2];
    nx = dx/L; ny = dy/L; nz = dz/L;
    #direction matrix
    I = [nx^2  nx*ny nx*nz;
         ny*nx ny^2  ny*nz;
         nz*nx nz*ny nz^2];

end
