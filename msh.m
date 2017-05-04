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
            tbe = [];%zeros(size(obj.Edges,1));
            for row = 1:size(obj.TetrahedronsByPoints, 1)
                tetrahedron = obj.TetrahedronsByPoints(row,:);
                pairs = combnk(tetrahedron,2);
                [~,indexes] = ismember(pairs, obj.Edges, 'rows');
                [~,tmp] = ismember(pairs, fliplr(obj.Edges), 'rows');
                indexes = indexes + tmp;
                tbe = [tbe; indexes'];
            end
        end
    end
    
end

