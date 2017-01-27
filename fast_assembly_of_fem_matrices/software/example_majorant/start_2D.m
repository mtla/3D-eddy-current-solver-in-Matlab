
fprintf('\n')
fprintf('Majorant calculation (2D) for a given approximation of the Poisson problem.\n')
fprintf('NOTE: Runtimes are reported in square brackets.\n')
fprintf('---------------------------------------------------------------------------\n\n')

clear all
close all
add_paths

levels     = 3;            % how many times the mesh is refined before solving
draw       = 1;            % 1=visualize; 0=no visualization
plots      = 2;            % 1=compact u_h,y_h ; 2=u_h, y_h, exact error, majorant indicator
testfolder = 'test2D';     % specify the test folder

% get Poisson problem data from the test folder
cd(testfolder)
load mesh;                 % load initial mesh
f = @loadfun;              % right hand side
u = @exact;                % exact solution
load const;                % load Cf, the Friedrichs constant
cd ..

for i=1:levels
    [nodes2coord,elems2nodes,bedges2nodes] = refinement_uniform_2D(nodes2coord,elems2nodes,bedges2nodes);
end

% calculate affine transformations and sign data for RT basis functions
%----------------------------------------------------------------------
tic;
[B_K,b_K,B_K_det] = affine_transformations(nodes2coord,elems2nodes);
elems2edges       = get_edges(elems2nodes);
signs             = signs_edges(elems2nodes);
time = toc;
fprintf('  Calculating affine transf., edge data, and edge signs: [%f]\n\n',time)

nelems = size(elems2nodes,1);
nnodes = size(nodes2coord,1);
nedges = max(max(elems2edges));

fprintf('  NOF elements  %d\n',nelems)
fprintf('  NOF nodes     %d\n',nnodes)
fprintf('  NOF edges     %d\n',nedges)

% solve the Poisson problem with zero Dirichlet BC
%-------------------------------------------------
[Su,Lu] = solver_poisson(elems2nodes,B_K,b_K,B_K_det,f);
tic;
bnodes2nodes    = unique(bedges2nodes);
dofs            = setdiff(1:nnodes,bnodes2nodes);
u_h             = zeros(nnodes,1);
u_h(dofs)       = Su(dofs,dofs) \ Lu(dofs);
time = toc;
fprintf('    Solving the linear system:  [%f]\n\n',time)

% calculate numerical error comparing with the exact solution
err = num_error( elems2nodes, B_K, b_K, B_K_det, u, u_h );
E = sqrt(sum(err));
fprintf('  Numerical error in energy norm: %f\n\n',E)

% calculate the majorant
%-----------------------
[Sy,My,L1y,L2y] = solver_majorant(elems2nodes, elems2edges, B_K, b_K, B_K_det, signs, f, u_h);

tic;
% these are needed for calculating the majorant value with matrix operations
uhSuh = u_h'*Su*u_h;
f_L2norm = norm_L2(elems2nodes,B_K,b_K,B_K_det,f);
beta = 1;
for i=1:10
    % solve the matrix system
    y_h = (  (1+beta) * Cf^2 * Sy  + (1+1/beta) * My  ) \ ...
          ( -(1+beta) * Cf^2 * L1y + (1+1/beta) * L2y );
    y_h = full(y_h);
          
	% calculate the two terms of the majorant with matrix operations
    %maj_residual = f_L2norm + y_h'*Sy*y_h + 2*y_h'*L1y;
    %maj_dual     = uhSuh + y_h'*My*y_h - 2*y_h'*L2y;
    
    % calculate the two terms of the majorant by integration
    [maj_residual,maj_dual] = majorant(elems2nodes, elems2edges, B_K, b_K, B_K_det, signs, f, u_h, y_h);
    
    % the majorant value
    maj = sqrt( (1+beta)*Cf^2*sum(maj_residual) + (1+1/beta)*sum(maj_dual) );
    
    % efficiency index of the majorant
    Ieff = maj/E;
    
    % print
    fprintf('    Iteration %d ',i);
    fprintf(': beta=%f ' ,beta);
    fprintf('; maj=%f '  ,maj );
    fprintf('; Ieff=%f\n',Ieff);
    
    % stopping criterion
    if ( i>1 && (maj_prev-maj)/maj_prev < 10^(-4) )
        break;
    end
    
    beta = sqrt(sum(maj_dual))/(Cf*sqrt(sum(maj_residual)));
    maj_prev = maj;
end
time = toc;
fprintf('    Time taken for iterations:  [%f]\n\n',time)


% visualize functions
%--------------------

if ( draw == 1 )
    %vangle = [-70.5,30];
    %vangle = 2;
    vangle = [-58,26];
    
    % approximation u_h
    %------------------
    [u_h_val,u_h_gval] = solution_P1(elems2nodes,B_K,u_h,1);
    u_h_gval_x = u_h_gval(:,1,:);
    u_h_gval_y = u_h_gval(:,2,:);
    
    % approximation y_h
    % note that the approximation is elementwise linear but the visualization
    % is done by taking a constant value (the midpoint) on each element
    %------------------------------------------------------------------------
    [y_h_val,y_h_dval] = solution_RT0(elems2edges,B_K,B_K_det,signs,y_h,1);
    y_h_val_x = y_h_val(:,1,:);
    y_h_val_y = y_h_val(:,2,:);
    
    if ( plots == 1 )
        
        figure
        show_nodal_scalar(u_h,nodes2coord,elems2nodes)
        axis on
        title('nodal approximation v of the Poisson problem')
        view([-58,26]);
        
        figure
        show_constant_scalar(y_h_val_x,nodes2coord,elems2nodes)
        axis on
        title('edge aproximation of the flux \tau (x-component)')
        view([-34,32]);

        figure
        show_constant_scalar(y_h_val_y,nodes2coord,elems2nodes)
        title('edge aproximation of the flux \tau (y-component)')
        axis on
        view([-53,34]);
        
    elseif ( plots == 2 )

        % u_h
        %----
        figure
        subplot(1,3,1);
        show_nodal_scalar(u_h,nodes2coord,elems2nodes)
        title('approx. u_h','FontWeight','bold')
        shading flat
        colorbar('Southoutside');
        axis off
        view(vangle);

        subplot(1,3,2);
        show_constant_scalar(u_h_gval_x,nodes2coord,elems2nodes)
        title('approx. grad(u_h) x-comp.','FontWeight','bold')
        shading flat
        colorbar('Southoutside');
        axis off
        view(vangle);

        subplot(1,3,3);
        show_constant_scalar(u_h_gval_y,nodes2coord,elems2nodes)
        title('approx. grad(u_h) y-comp.','FontWeight','bold')
        shading flat
        colorbar('Southoutside');
        axis off
        view(vangle);

        % y_h
        %----
        figure
        subplot(1,3,1);
        show_constant_scalar(y_h_val_x,nodes2coord,elems2nodes)
        title('approx. y_h x-comp.','FontWeight','bold')
        shading flat
        colorbar('Southoutside');
        axis off
        view(vangle);

        subplot(1,3,2);
        show_constant_scalar(y_h_val_y,nodes2coord,elems2nodes)
        title('approx. y_h y-comp.','FontWeight','bold')
        shading flat
        colorbar('Southoutside');
        axis off
        view(vangle);

        subplot(1,3,3);
        show_constant_scalar(y_h_dval,nodes2coord,elems2nodes)
        title('approx. div(y_h)','FontWeight','bold')
        shading flat
        colorbar('Southoutside');
        axis off
        view(vangle);

        % error
        %------
        figure
        subplot(1,2,1)
        areas = 1/2.*abs(B_K_det);
        show_constant_scalar(err./areas,nodes2coord,elems2nodes)
        title('exact error distribution')

        subplot(1,2,2)
        show_constant_scalar(((1+1/beta).*maj_dual)./areas,nodes2coord,elems2nodes)
        title('distribution of dual part of majorant')
        
    end
    
end
