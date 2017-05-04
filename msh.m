classdef msh
    %MSH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        DT
        Points
        ConnectivityList
        Edges
        tetrahedron_edges
    end
    
    methods
        function obj = msh(vertices)
            obj.DT = delaunayTriangulation(vertices);
            obj.Points = obj.DT.Points;
            obj.ConnectivityList = obj.DT.ConnectivityList;
            obj.Edges = edges(obj.DT);
        end
        function tetramesh(obj)
            tetramesh(obj.DT)
        end
    end
    
end

