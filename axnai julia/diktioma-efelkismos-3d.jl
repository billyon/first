using Plots
#komvos1 komvos2 komvos3
xyz = 0.1*[0 1 -1;
           0 2 -1;
           0 0  0
]; #orismos komvon san xyz
n = size(xyz)[2]; #arithmos komvon

#komvos1_start komvos1_end
numel = [1 2;
         1 3
]; #orismos ton stixion os grammes me toys akraious komvous

#orismos elastikotitas stixion
E = 220e9.*ones(N,1); A = 4e-4.*ones(N,1);

N = size(numel)[1]; #arithmos stoixion


#orismos ton desmeumenon komvon
BC = zeros(3*n,1)
#for 2d uncomment below line
dimconstr = 1:3:3*n; BC[dimconstr] .= 1;
BC[4:end] .= 1; #0==free
fixedofs = findall(x->x==1,BC[:])
freedofs = findall(x->x==0,BC[:])
#boundary conditions
U = zeros(3*n,1); U[4:end].=0;
U[dimconstr].=0; #uncomment for 2d
F = zeros(3*n,1); F[[1 2]]=[-100 100] ;
#orismos olikou mitroou stivarotitas
global K = zeros(3*n,3*n);
#start main loop to calculate stiffness matrix
for i=1:N
    i1 = numel[i,1]
    i2 = numel[i,2]
    #calclulate directions
    x1 = xyz[1,i1]; x2 = xyz[1,i2]; dx = x2-x1;
    y1 = xyz[2,i1]; y2 = xyz[2,i2]; dy = y2-y1;
    z1 = xyz[3,i1]; z2 = xyz[3,i2]; dz = z2-z1;
    L = sqrt(dx^2+dy^2+dz^2);
    nx = dx/L; ny = dy/L; nz = dz/L;
    #direction matrix
    I = [nx^2  nx*ny nx*nz;
         ny*nx ny^2  ny*nz;
         nz*nx nz*ny nz^2];
    k = E[i]*A[i]/L
    Ke = k*[I -I; -I I]
    #println(size(Ke),"\n")
    edofs = [3*i1-2 3*i1-1 3*i1  3*i2-2 3*i2-1 3*i2]
    K[edofs,edofs] = K[edofs,edofs][1,:,1,:] + Ke;
    #println(K,"\n")
end
Kfinal = K[freedofs,freedofs]
rhs = F
rhs = rhs[freedofs]
rhs = rhs - K[freedofs,fixedofs]*U[fixedofs]

Ufree = Kfinal\rhs
println(Ufree)
