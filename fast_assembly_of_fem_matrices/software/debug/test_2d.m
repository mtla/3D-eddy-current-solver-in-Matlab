%####################################################################
% This function tests calc_meshdata by integrating a function over
% the unit square in 2D. The domain is divided into triangles, and
% the integration quadrature is defined on the unit triangle, so the
% integration over the domain is performed by chage of variables.
% NOTE: This function performs the integration both by a loop over
% elements (the 'normal' way), and in a vectorized way, which is
% much faster.
%####################################################################

function test_2d
    add_paths

    % handle mesh and affine mappings
    cd test2D
    load mesh;                 % load initial mesh
    cd ..
    for i=1:7
        [nodes2coord,elems2nodes,bedges2nodes] = refinement_uniform_2D(nodes2coord,elems2nodes,bedges2nodes);
    end
    nelems = size(elems2nodes,1);
    disp(['Number of elements: ',num2str(nelems)])
    [B_K,b_K,B_K_det] = affine_transformations(nodes2coord,elems2nodes);

    % integration quadrature in 2d
    [ip,w,nip] = intquad(4,2);
    
    % calculate the integral of the function 'fun' in a normal way,
    % i.e., by a loop over elements
    tic
    integral = 0;
    for i=1:nelems
        % get element info
        B_K_     = B_K(:,:,i);
        b_K_     = b_K(i,:)';
        B_K_detA = abs(B_K_det(i));

        % calculate the affine mapping F_K( ip )
        F_K_ip = B_K_ * ip' + b_K_(:,ones(1,nip));
        
        % the values of the function to be integrated
        funval = fun(F_K_ip');
        
        % calculate the integral contribution
        integral = integral + sum( w .* B_K_detA .* funval );
    end
    time = toc;
    disp(['Time taken to calculate integral (normal):     [',num2str(time),']'])
    disp(['The integral value (normal):     ',num2str(integral)])
    
    
    % calculate the integral of the function 'fun' in a vectorized way,
    % i.e., by a loop over integration points
    tic;
    integral = zeros(nelems,1);
    B_K_detA = abs(B_K_det);
    for i=1:nip
        % calculate the affine mapping F_K( ip )
        F_K_ip = squeeze(amsv(B_K, ip(i,:)))' + b_K;

        % the values of the function to be integrated
        funval = fun(F_K_ip);

        % calculate the integral contribution
        integral = integral + w(i) .* B_K_detA .* funval;
    end
    integral = sum(integral);
    time = toc;
    disp(['Time taken to calculate integral (vectorized): [',num2str(time),']'])
    disp(['The integral value (vectorized): ',num2str(integral)])
    
end


%####################################################################
% this is the function we are integrating
%####################################################################

function val = fun(p)
    %val = ones(size(p,1),1);    %fun = 1
    %val = p(:,1);               %fun = x
    val = p(:,2).^2;            %fun = y^2
end