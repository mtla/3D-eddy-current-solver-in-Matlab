function test_average()

add_paths
cd test2D
load mesh;                 % load initial coarse mesh

for i=1:2  %refine mesh uniformly
    [nodes2coord,elems2nodes] = refinement_uniform_2D(nodes2coord,elems2nodes);
end

cd ..
[B_K,b_K,B_K_det] = affine_transformations(nodes2coord,elems2nodes);
B_K_detA = abs(B_K_det);
f = @fun;

disp(['Number of elements: ',num2str(size(elems2nodes,1))])

avgs = average(elems2nodes,B_K,b_K,B_K_det,f);  %computing averages using the default quadrature rule
integral = sum(avgs.*(B_K_detA./2));
fprintf('integral of the function computed from the element averages = %18.16f (default quadrature)\n',integral)

for qr=3:12 %quadrature rule
    avgs = average(elems2nodes,B_K,b_K,B_K_det,f,qr);  %computing averages using the default quadrature rule
    integral = sum(avgs.*(B_K_detA./2));
    fprintf('integral of the function computed from the element averages = %18.16f (order of quadrature = %d)\n',integral,qr)
end

%figure(1); show_mesh(elems2nodes,nodes2coord)
figure(2); show_constant_scalar(avgs,nodes2coord,elems2nodes)

end


% the function whose average we want to calculate
function val = fun(p)
    val = exp( p(:,1).^2 + p(:,2).^2 );
    %val = ones(size(p,1),1);    %fun = 1
    %val = p(:,1);               %fun = x
    %val = p(:,2).^2;            %fun = y^2
    
    %the exact value of integral over the unit box is
    %2.139350129805326552121276528614617425580337155525191298143120...
end