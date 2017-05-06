classdef msh
    % We just create our own class and add functionality to it because 
    % matlab does not allow us to extend the delaunayTriangulation
    % class.
    
    properties
        DT
        Points
        TetrahedronsByPoints
        Edges
        TetrahedronsByEdges
    end
    
    methods
        function obj = msh(vertices)
            obj.DT = delaunayTriangulation(vertices);
            obj.Points = obj.DT.Points;
            obj.TetrahedronsByPoints = obj.DT.ConnectivityList;
%             obj.TetrahedronsByPoints = obj.removeDuplicateTetrahedrons();
            obj.Edges = edges(obj.DT);
            obj.TetrahedronsByEdges = tetrahedrons2edges(obj);
        end
        function tetramesh(obj)
            tetramesh(obj.TetrahedronsByPoints, obj.Points)
%             tetramesh(obj.DT)
        end
        function plot3(obj, values)
           % normalize values
           if (exist('values','var')==1)
               range = max(values) - min(values);
               values = (values - min(values)) / range;
           end
           for ne = 1:obj.ne()
               hold on
               e = obj.getEdgeCoordinates(ne);
               if (exist('values','var')==1)
                   color = [values(ne) 0 1-values(ne)];
                   plot3(e(:,1),e(:,2),e(:,3),'Color', color);
               else
                   plot3(e(:,1),e(:,2),e(:,3));
               end
           end
        end
        function scatter3(obj, values)
           % normalize values
           if (exist('values','var')==1)
               range = max(values) - min(values);
               values = (values - min(values)) / range;
           end
           for ne = 1:obj.np()
               hold on
               node = obj.Points;
               if (exist('values','var')==1)
                   color = [values(ne) 0 1-values(ne)];
                   scatter3(node(:,1),node(:,2),node(:,3),'MarkerFaceColor', color);
               else
                   scatter3(node(:,1),node(:,2),node(:,3));
               end
           end
        end
        function p = points(obj)
            p = obj.Points;
        end
        function e = edges(obj)
            e = obj.Edges;
        end
        function np = np(obj) % number of points in mesh
            np = size(obj.Points,1);
        end
        function ne = ne(obj) % number of edges in mesh
            ne = size(obj.Edges,1);
        end
        function nt = nt(obj) % number of tetrahedrons in mesh
            nt = size(obj.TetrahedronsByPoints, 1);
        end
        function edge = getEdge(obj, edgeID)
            edge = obj.Edges(edgeID,:);
        end
        function coordinates = getEdgeCoordinates(obj, edgeID)
            edge = obj.getEdge(edgeID);
            coordinates = obj.Points(edge,:);
        end
        function points = tetrahedron2points(obj, tetrahedronID)
            points = obj.TetrahedronsByPoints(tetrahedronID,:);
        end
        function edges = tetrahedron2edges(obj, tetrahedronID)
            tetrahedron = obj.TetrahedronsByPoints(tetrahedronID,:);
            pairs = combnk(tetrahedron,2);
            [~,b] = ismember(pairs, obj.Edges, 'rows');
            [~,tmp] = ismember(pairs, fliplr(obj.Edges), 'rows');
            edges = b' + tmp';
        end
%     end
%     methods(Access = private)

        function tbe = tetrahedrons2edges(obj)
            n = obj.nt(); % number of tetrahedrons
            tbe = zeros(n, 6); 
            for row = 1:n
                tbe(row,:) = obj.tetrahedron2edges(row);
            end
        end
        function tetrahedrons = removeDuplicateTetrahedrons(obj)
            tetrahedrons = obj.DT.ConnectivityList;
%             a = (1:obj.np())';
%             occurence = histc(obj.TetrahedronsByPoints(:),a);
            for n = 1:size(tetrahedrons,1)
%                 n
%                 obj.ne()
                [rows,~] = find(tetrahedrons==n,size(tetrahedrons,1),'first');
                tetrahedrons(rows(5:end),:) = [];
            end
        end
    end
    
end

