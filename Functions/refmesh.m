%refmesh refine the mesh; update the current density and reluctivity.
% 
% Borrowed and adapted from the mathematics department course on FEM plus
% the SMEKlib opensource library.
% Copyright Antti Hannukainen 2007, Antti Lehikoinen 2013-2017


function [rmesh, currentDensity, reluctivity] = refmesh(msh, currentDensity, reluctivity )

msh = inittri(msh.p, msh.t);

t = msh.t;
p = msh.p;
e = msh.edges;
t2e = msh.t2e;

Nt = size(t,2);
Ne = size(e,2);

for n=1:2
  e_nodes(n,:) = sum( [p(n,e(1,:)) ; p(n,e(2,:))])/2;
end

% Create new mesh
rmesh.p = [p  e_nodes];


% Create New elements
eb = size(p,2); 

% Edges as n1->n2, n2->n3, n1->n3.
new_t = [t(1,:) ; t2e(1,:)+eb ; t2e(3,:)+eb ];
rmesh.t = [new_t];

new_t = [t(2,:) ; t2e(1,:)+eb ; t2e(2,:)+eb ];
rmesh.t = [rmesh.t  new_t];

new_t = [t(3,:) ; t2e(3,:)+eb ; t2e(2,:)+eb ];
rmesh.t = [rmesh.t  new_t];

new_t = [t2e(1,:)+eb ; t2e(2,:)+eb ; t2e(3,:)+eb ];
rmesh.t = [rmesh.t  new_t];

rmesh = struct('t', rmesh.t, 'p', rmesh.p);

% updating element properties
currentDensity = repmat(currentDensity, 1, 4);
reluctivity = repmat(reluctivity, 1, 4);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following functions taken from the Open Source SMEKlib library
% (c) 2013-2017 Antti Lehikoinen / Aalto University

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function msh = inittri(p, t)
%inittri initializes a mesh struct
% 
% msh = inittri(p, t) initializes a mesh from nodes p and elements t
% The struct contains the following fields
%   p = nodes as 2xNp array
%   t = elements as 3xNe array
%   edges = edges as 2xNedges array
%   t2e = 3xNe array describing which edges each element has
%   e2t = 2xNedges array describing which elements each edge belongs to
%
% Copyright (c) 2016 Antti Lehikoinen / Aalto University,
%   based on the work of Antti Hannukainen and Mika Juntunen, 3.11.2005

t = sort(t, 1);

%defining edges an' pals
[edges, e2t, t2e] = getEdges(t);

%struct
msh = struct('t', t, 'p', p, 'edges', edges, 'e2t', e2t, 't2e', t2e);

end

function [edges, e2t, t2e] = getEdges(t)
%getEdges returns the edge definition.
%
% [edges, e2t, t2e] = getEdges(t)
% when given a 3xNe triangulation t, the function returns the following
% arrays
%   edges = 2xNedges array, containing the start- and end-nodes of each
%      edge in the mesh (on first and second row respectively)
%   e2t = a 2xNedges array, containing the triangles that the mesh belongs
%       to
%   t2e = a 3xNe array, containing the edges that each triangle has
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University,
%   based on the work of Antti Hannukainen and Mika Juntunen, 3.11.2005

%

% initialization
t = sort(t,1);
Nt = size(t, 2);

%edges
edges = sort([t(1,:) t(2,:) t(3,:);
    t(2,:) t(3,:) t(1,:)], 1);

[edges, ~, t2e] = sc_unique(edges);
t2e = reshape(t2e,Nt,3)';

%getting the edges-to-elements array
e = [t2e(1,:) t2e(2,:)  t2e(3,:)];
t_list = repmat(1:Nt, 1,3);

[ef,If]= unique(e, 'first');
[el,Il]= unique(e, 'last');

e2t(1,ef) = t_list(If);
e2t(2,el) = t_list(Il);

e2t(2, (e2t(1,:)-e2t(2,:))==0)=0;

end

function [B,I,J] = sc_unique(A)

% sort columnwise
A = sort(A, 1);

% unique
[B,I,J] = unique(A','rows');

%untranspose
B = B';

%unrolling
I = I(:)';
J = J(:)';

end