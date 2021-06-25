function kspace = parse_list_file(filename)
% PARSE_LIST_FILE parses Philips list file without touching the .data file
% can deal with spaces in property names (as seen in R5.3.31)
% based on loadRawKspacePhilips.m by Wouter Potters (w.v.potters@amc.nl)
%
% Andreas Wetscherek (a.wetscherek@icr.ac.uk), 06/07/2020

fid = fopen([filename(1:end-4) 'list']);
currentline = '#';
kspace_properties = struct();
while any(strcmp(currentline(1),{'#','.'}))
    currentline = fgetl(fid);
    if strcmp(currentline(1),'.')
        % save property;
        C = textscan(currentline,'. %f %f %f %s : %f %f',1);
        index = 4;
        if ~isempty(C{5})
            fieldname = strrep(C{4}{1},'-','_');
        else
            C = textscan(currentline,'. %f %f %f %s %s %s %s : %f %f',1);
            fieldname = strrep(C{4}{1},'-','_');
            index = 4;
            while index < 7 && ~strcmp(C{index + 1}{1}, ':')
                index = index + 1;
                fieldname = [fieldname '_' strrep(C{index}{1},'-','_')];
            end
            if index < 7 && strcmp(C{index + 1}{1}, ':')
                C = textscan(currentline,'. %f %f %f %s %s %s : %f %f',1);
            end
        end
        kspace_properties.(fieldname) = [C{index+1} C{index+2}];
    end
end

fseek(fid, -numel(currentline), 'cof');

list = textscan(fid,'%s %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',inf,'MultipleDelimsAsOne',1);
fclose(fid);

while strcmp(list{:,1}(end),'#')
    for i = length(list):-1:1
        try
            list{:,i}(length(list{:,1})) = [];
        catch
            %ignore error
        end
    end
end

headers = 'typ mix   dyn   card  echo  loca  chan  extr1 extr2 ky    kz    n.a.  aver  sign  rf    grad  enc   rtop  rr    size   offset';
headers = textscan(headers,'%s');
headers = headers{1}(:)';

for i = 1:length(headers)
    eval(['kspace.' strrep(headers{i},'.','_') '= list{:,' num2str(i) '};'])
end

kspace.kspace_properties = kspace_properties;