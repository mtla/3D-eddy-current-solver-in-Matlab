classdef msh < handle
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
%             tetrarefine3(obj);
%             obj.TetrahedronsByPoints = obj.removeDuplicateTetrahedrons();
            obj.Edges = edges(obj.DT);
            obj.TetrahedronsByEdges = tetrahedrons2edges(obj);
        end
        function tetramesh(obj)
            tetramesh(obj.TetrahedronsByPoints, obj.Points)
        end
        function tetrarefine3(obj)
            [p, t,~] = tetrarefine3(obj.Points, obj.TetrahedronsByPoints);
            obj.Points = p;
            obj.TetrahedronsByPoints = t;
        end
        function plot3(obj, values)
           % normalize values
           if (exist('values','var')==1)
               range = max(values) - min(values);
               values = (values - min(values)) / range - min(values);
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
           hold on
           node = obj.Points;
           % normalize values
           if (exist('values','var')==1)
               range = max(values) - min(values);
               values = (values - min(values)) / range - min(values);
               color = [values zeros(obj.np(),1) 1-values];
               
               scatter3(node(:,1),node(:,2),node(:,3), 50, color);%color);
           else
               scatter3(node(:,1),node(:,2),node(:,3));
           end
        end
        function p = points(obj)
            p = obj.Points;
        end
        function [e, all] = edges(obj)
            n = obj.nt(); % number of tetrahedrons
            e = zeros(n*6,2);
            for tID = 1:n
                t = obj.TetrahedronsByPoints(tID,:); % tetrahedron
                edges = [t(2) t(4); ...
                         t(2) t(3); ...
                         t(2) t(1); ...
                         t(4) t(3); ...
                         t(3) t(1); ...
                         t(1) t(4)];
                e(tID*6-5:tID*6,:) = edges;
            end
%             e
            all = e;
            e = unique(e,'rows')
            for en = 1:size(e,1)
                try
                    edge = e(en,:)
                catch
                end
%                 ismember(edge,e,'rows')
                [rows,~] = find(e==edge,1)
%                 [~, rows] = ismember(edge,e,'rows')
                e(rows(2:end),:) = [];
                [rows,~] = find(e==fliplr(edge),1)
                e(rows(2:end),:) = [];
            end
            
            commonRows = ismember(e,fliplr(e),'rows')
            e
            size(e)
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
%             t = obj.TetrahedronsByPoints(tetrahedronID,:);
%             edges = [t(2) t(4); ...
%                      t(2) t(3); ...
%                      t(2) t(1); ...
%                      t(4) t(3); ...
%                      t(3) t(1); ...
%                      t(1) t(4)];
        end
%     end
%     methods(Access = private)

        function tbe = tetrahedrons2edges(obj)
            [edges, tedges] = obj.edges();
            tbe = zeros(obj.nt(), 6); 
            for tID = 1:obj.nt()
                current_tetrahedron = tedges(tID*6-5:tID*6,:);
                tetrahedron = zeros(1,6);
                for i = 1:6
                    cedge = current_tetrahedron(i,:);
                    [~, tetrahedron(i)] = ismember(cedge,edges,'rows');
                end
                tbe(tID,:) = tetrahedron;
            end
        end
        function tetrahedrons = removeDuplicateTetrahedrons(obj)
            tetrahedrons = obj.DT.ConnectivityList;
%             a = (1:obj.np())';
%             occurence = histc(obj.TetrahedronsByPoints(:),a);
            for n = 1:size(tetrahedrons,1)
                [rows,~] = find(tetrahedrons==n,size(tetrahedrons,1),'first');
                tetrahedrons(rows(5:end),:) = [];
            end
        end
    end
    
end

