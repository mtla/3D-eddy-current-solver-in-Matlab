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
        end
        function tetramesh(obj)
            tetramesh(obj.DT)
        end
    end
    
end

