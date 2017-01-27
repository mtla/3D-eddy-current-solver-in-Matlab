%####################################################################
% This function tests calc_meshdata by integrating a function over
% the unit square in 3D. The domain is divided into tetras, and
% the integration quadrature is defined on the unit tetra, so the
% integration over the domain is performed by chage of variables.
% NOTE: This function performs the integration both by a loop over
% elements (the 'normal' way), and in a vectorized way, which is
% much faster.
%####################################################################

function test_3d
    add_paths

    [p,f,t] = read_tetgenmesh('3d/unitcube.5');
    nelems = size(t,2);
    disp(['Number of elements: ',num2str(nelems)])

    tic
    md = calc_meshdata(3,p,f,t);
    time = toc;
    disp(['Time taken to calculate meshdata:              [',num2str(time),']'])

    % integration quadrature in 3d
    [ip,w] = inttet(2);
    nip = size(ip,2);

    % calculate the integral of the function 'fun'
    tic
    integral = 0;
    for i=1:nelems
        % get element info
        B_K      = md.B_K(:,:,i);
        b_K      = md.b_K(:,i);
        B_K_detA = abs(md.B_K_det(i));

        % calculate the affine mapping F_K( ip )
        F_K_ip = B_K * ip + b_K(:,ones(1,nip));
        
        % the values of the function to be integrated
        funval = fun(F_K_ip);
        
        % calculate the integral contribution
        integral = integral + sum( w' .* B_K_detA .* funval );
    end
    time = toc;
    disp(['Time taken to calculate integral (normal):     [',num2str(time),']'])
    disp(['The integral value (normal):     ',num2str(integral)])
    
    
    % calculate the integral of the function 'fun' in a vectorized way,
    % i.e., by a loop over integration points
    tic;
    integral = zeros(1,nelems);
    B_K_detA = abs(md.B_K_det);
    for i=1:nip
        % calculate the affine mapping F_K( ip )
        F_K_ip = squeeze(amsv(md.B_K, ip(:,i))) + md.b_K;

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
% This is the function we are integrating.
%####################################################################

function val = fun(p)
    %val = ones(1,size(p,2));    %fun = 1
    %val = p(1,:);               %fun = x
    val = p(2,:).^2;            %fun = y^2
end