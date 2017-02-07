%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Describing the problem in EXAMPLE 1:
% defining the mesh, boundary conditions, source terms (current density)
% and material characteristics

%coordinates of the nodes in a 2x4 array [m]
p = [0 0;
    1 0;
    1 1;
    0 1]';

%elements as a 3x2 array
t = [1 2 3;
    1 3 4]';

%mesh struct that will be used by the provided matrix assembly functions
msh = struct('p', p, 't', t);

%determining nodes on the Dirichlet boundary
DirichletNodes = [2 3];%find( msh.p(1,:) > 0.99 );
DirichletNodes = unique(DirichletNodes); %ALWAYS make sure you don't have dublicates; they frak up the entire problem

%current density in each element [A/m^2]
currentDensity = [0 1];

%reluctivity of each element [A/(Tm)]
reluctivity = [1 1] / (pi*4e-7);

% example on refining the mesh
%[msh, currentDensity, reluctivity] = refmesh(msh, currentDensity, reluctivity);
%DirichletNodes = find( msh.p(1,:) > 0.99 );
%DirichletNodes = unique(DirichletNodes); %ALWAYS make sure you don't have dublicates; they frak up the entire problem


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solving problem
Np = size(msh.p,2); %number of nodes in the problem
freeNodes = setdiff(1:Np, DirichletNodes); %nodes NOT in Dirichlet bnd

%assembling stiffness matrix and load vector
[S,f] = calculate_Arrays(msh, reluctivity, currentDensity);

%calculating potentials in the free nodes
Afree = S(freeNodes,freeNodes) \ f(freeNodes);
% NOTE: this is equivalent to Afree = inv(S) * f, but much faster

%assembling solution in the entire region
A_total = zeros(Np,1);
A_total(freeNodes) = Afree;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postprocessing

%plotting mesh
figure(1); clf; hold on;
triplot(msh.t', msh.p(1,:), msh.p(2,:))
plot(msh.p(1, DirichletNodes), msh.p(2, DirichletNodes), 'ko');
xlabel('x [m]')
ylabel('y [m]')
title('The mesh used; Dirichlet nodes')

%plotting non-zero entries of S
figure(2); clf;
spy(S);
xlabel('Column, j')
ylabel('Row, i')
title('Non-zero entries of the stiffness matrix')

%plotting values of the vector potential A
figure(3); clf;
trimesh(msh.t', msh.p(1,:), msh.p(2,:), A_total);
xlabel('x [m]')
ylabel('y [m]')
title('Vector potential')

%drawing flux lines = equipotential lines of A plus the flux density
figure(4); clf; box on;
plotFluxDensity(A_total, msh);
draw_Equipotentials(A_total, msh, 5)

xlabel('x [m]')
ylabel('y [m]')
title('Flux lines and flux density')