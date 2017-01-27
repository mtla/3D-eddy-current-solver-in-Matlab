function show_nodal_scalar(nodalValue,nodes2coord,elems2nodes,nodalDisplacement)
switch nargin, 
    case 3,
        nodalDisplacement=0*nodes2coord;
    case {0, 1, 2}
        fprintf('missing parameters')
end


X=nodes2coord(:,1)+nodalDisplacement(:,1);
Y=nodes2coord(:,2)+nodalDisplacement(:,2);

%fill3(X(elems2nodes)',Y(elems2nodes)',nodalValue(elems2nodes)',nodalValue(elems2nodes)');
fill3(X(elems2nodes)',Y(elems2nodes)',nodalValue(elems2nodes)',nodalValue(elems2nodes)','FaceColor','interp','LineStyle','none');
set(gcf,'renderer','zbuffer');