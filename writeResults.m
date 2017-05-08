function writeResults( results )
%WRITERESULTS Summary of this function goes here
%   Detailed explanation goes here

    [filename, filefolder] = uiputfile({'*.csv','Comma separated file';'*.txt','Raw Text File';'*.vtk','Binary Paraview file'}, 'Save output');
    output_path = strcat(filefolder,filename);
    
    [~, ~, fileformat] = fileparts(output_path);
    
    switch fileformat
        case '.csv'
            fid = fopen(output_path, 'w');
            fprintf(fid, 'X coords,Y coords,Z coords,Scalar\n');
            fclose(fid);
            
            dlmwrite(output_path, results, '-append', 'precision', '%.6f', 'delimiter', ',');
%             csvwrite(output_path, results)
        otherwise
            disp(strcat('Output format  ', fileformat, ' is not supported (yet?)'));
    end
end

