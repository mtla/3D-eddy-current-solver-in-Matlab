function [ global_matrix ] = element2global( tetrahedron_matrix, N_matrix)
% ELEMENT2GLOBALMATRIX This function takes a 4x4 tetrahedron matrix and nxn matrix 
% and assembles a global matrix according to the nxn matrix
% 
% Input:
%     tetrahedron_matrix:
%         4x4 matrix with values according to the shape function
%     N_matrix:
%         "a recepi" to add the tetrahedron_matrix elements to global_matrix positions
%         N_matrix's assumed columns: 1) element number 2) number of nodes in element (4) 
%         3) 1st global numbers of nodes 4) 2nd global numbers of nodes ... 6) 4th global numbesr of nodes
%         7) material number 8) source number
% 
% Output:
%     global_matrix:
%         Global matrix with values from tetrahedron_matrix according to he N_matrix
% Syntax: element2global(tetrahedron_matrix, N_matrix)
%

% global_matrix is first created as a matrix of zeros 
% size depends on the largest value of N_matrix
% NOT SURE YET IF CORRECT!!
global_matrix = zeros(max(max(N_matrix)),max(max(N_matrix)));

for n=1:7
    for x=1:4
         for y=1:4
              % 
              global_matrix(N_matrix(n, x+2), N_matrix(n, y+2)) = tetrahedron_matrix(x,y);
         end
    end
end
