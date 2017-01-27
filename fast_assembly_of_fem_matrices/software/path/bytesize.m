function str = bytesize(in, fid)

% BYTESIZE
%	Writes the memory usage of the provide variable to the given
%	file identifier. Output is written to screen if fid is 1, empty
%   or not provided.
%
%   Code by "MatlabSorter" in the forum StackOverflow.com :
%   http://stackoverflow.com/questions/4845561/how-to-know-the-size-of-a-variable-in-matlab
%	(17.09.2014)
%

if nargin == 1 || isempty(fid)
    fid = 1;
end

s = whos('in');
str = Bytes2str(s.bytes);
%fprintf(fid,[str,'\n']);

end


%#####################################################################
%#####################################################################


function str = Bytes2str(NumBytes)
% BYTES2STR
%	Private function to take integer bytes and convert it to
% 	scale-appropriate size.

scale = floor(log(NumBytes)/log(1024));
switch scale
    case 0
        str = [sprintf('%.0f',NumBytes) ' b'];
    case 1
        str = [sprintf('%.2f',NumBytes/(1024)) ' kb'];
    case 2
        str = [sprintf('%.2f',NumBytes/(1024^2)) ' Mb'];
    case 3
        str = [sprintf('%.2f',NumBytes/(1024^3)) ' Gb'];
    case 4
        str = [sprintf('%.2f',NumBytes/(1024^4)) ' Tb'];
    case -inf
        % Size occasionally returned as zero (eg some Java objects).
        str = 'Not Available';
    otherwise
       str = 'Over a petabyte!!!';
end
end