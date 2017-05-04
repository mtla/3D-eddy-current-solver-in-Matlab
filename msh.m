classdef msh
    % Because matlab does not allow us to extend the delaunayTriangulation
    % class we just create our own class and add functionality to it.
    
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
            obj.Edges = edges(obj.DT);
            obj.TetrahedronsByEdges = tetrahedrons2edges(obj);
        end
        function tetramesh(obj)
            tetramesh(obj.DT)
        end
        function p = points(obj)
            p = obj.Points;
        end
        function e = edges(obj)
            e = obj.Edges;
        end
%         function s = size(obj)
%             s = [size(obj.Points); size(obj.Edges); size(obj.TetrahedronsByPoints)];
%         end
    end
    methods(Access = private)
        function tbe = tetrahedrons2edges(obj)
            n = size(obj.TetrahedronsByPoints, 1); % number of tetrahedrons
            tbe = zeros(n, 6); % ugly. What could the number of edges 
            for row = 1:n
                tetrahedron = obj.TetrahedronsByPoints(row,:);
                pairs = combnk(tetrahedron,2);
                [~,b] = ismember(pairs, obj.Edges, 'rows')
                [~,tmp] = ismember(pairs, fliplr(obj.Edges), 'rows')
                tbe(row,:) = b' + tmp';
            end
        end
    end
    
end

