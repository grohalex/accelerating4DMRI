function kspace = loadRawKspacePhilips(listordatafile)
% LOADRAWKSPACEPHILIPS Loads the Philips *.list and *.data files that can
% be exported from Philips MRI systems.
%
% [KSPACE] = LOADRAWKSPACEPHILIPS(LISTORDATAFILE)
%
% KSPACE          Struct containing all kspace parameters read from the 
%                 header file (.list) and binary data file (.data).
%                 More information see below.
% LISTORDATAFILE  valid .list or .data file. Note that the data and list file should be 
%                 located in the same folder.
%
% kspace.typ          type kline (STD,REJ,NOI)
% kspace.mix          mixed sequence number
% kspace.dyn          dynamic scan number
% kspace.card         cardiac phase number
% kspace.echo         echo number
% kspace.loca         location number
% kspace.chan         synco channel number
% kspace.extr1        extra attribute 1
% kspace.extr2        extra attribute 2
% kspace.[ky,kz]      k-space location in 1st and 2nd preparation direction
%                      (image data)
% kspace.[kx,ky,kz]   k-space location in 1st, 2nd and 3rd preparation direction
%                      (spectroscopy data)
% kspace.n_a_         ?
% kspace.aver         sequence number of this signal average
% kspace.sign         sign of measurement gradient used for this data
%                      vector (1=pos, -1=neg)
% kspace.rf           sequence number of this rf echo
%                       (only for TSE,TFE,GraSE)
% kspace.grad         sequence number of this gradient echo
%                       (only for EPI, GraSE)
% kspace.enc          encoding time (only for EPI, GraSE)
% kspace.rtop         R-top offset in ms
% kspace.rr           RR interval in ms
% kspace.size         data vector size in bytes (1 complex elements = 2
%                       single floats = 8 bytes)
% kspace.offset       data vector offset in bytes (first data vector starts
%                       at offset 0)
% kspace.complexdata  calculated complex kspace data (kx complex numbers at location ky,kz)
%
% kspace.kspace_properties  Additional values found in the header file, such
%                            as kspace offsets.
%
% For more information on the header, please refer to the *.list file using an arbitrary text editor.
%
%
% Example:
% 
% kspace = loadRawKspacePhilips('C:/raw_004.list')
%
% Wouter Potters 
% Academic Medical Center
% Amsterdam, The Netherlands
% bug found? question? mail me at w.v.potters@amc.nl
% 

if ischar(listordatafile)
    switch listordatafile(end-4:end)
        case '.list'
            listfile = listordatafile;
            datafile = [listordatafile(1:end-5) '.data'];
        case '.data'
            datafile = listordatafile;
            listfile = [listordatafile(1:end-5) '.list'];
        otherwise
            error('MATLAB:loadRawKspace:wrong_inputfile','Please provide a valid list or data file');
    end
end

% load the listfile
kspace = loadListFile(listfile);
kspace = loadDataFile(datafile,kspace);



function kspace = loadDataFile(datafile,kspace)
% kspace = loadDataFile(datafile,kspace)
kspace.complexdata = []; %create extra variable.
nroflocations = numel(kspace.typ);

% read NOI_binary_data
try 
    fid = fopen(datafile);
catch %#ok<*CTCH>
    error('MATLAB:loadDataFile:invalidFile','Data file is invalid.');
end
for cur = 1:nroflocations
    kspace.complexdata{cur,1} = [1 1i] * fread(fid,[2 kspace.size(cur)/8],'single',0,'ieee-le');
    %if ~rem(cur,10000)
    %    disp([num2str(100*(double(kspace.offset(cur))/sum(kspace.size))) '%'])  %8 = 4bytes*2rows(im+re)
    %end
end

fgetl(fid); %go one line further to reach end-of-file

if feof(fid)
    fclose(fid);
else
    n = 0;
    fprintf('=============================================\n')
    while ~feof(fid)
        fgetl(fid)
        n=n+1;
    end
    warning('MATLAB:loadDataFile:DidNotReachEnd',['End of file finally reached after ' num2str(n) ' extra lines. Please make sure all data is available.']);
    fprintf('=============================================\n')
end

function kspace = loadListFile(listfile)
%kspace = loadListfile(filepath)
headers = 'typ mix   dyn   card  echo  loca  chan  extr1 extr2 ky    kz    n.a.  aver  sign  rf    grad  enc   rtop  rr    size   offset';
headers = textscan(headers,'%s');
headers = headers{1}(:)';
[nrofheaderlines, kspace_properties] = getNrofheaderlines(listfile);

try 
    fid = fopen(listfile);
catch
    error('MATLAB:loadListFile:invalidFile','List file is invalid.');
end


list = textscan(fid,'%s %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',inf,'HeaderLines',nrofheaderlines,'MultipleDelimsAsOne',1);
fclose(fid);

% remove
while strcmp(list{:,1}(end),'#')
    for i = length(list):-1:1
        try
            list{:,i}(length(list{:,1})) = [];
        catch
            %ignore error
        end
    end
end

kspace = struct(); %create empty struct
%fill struct with information on klines
for i = 1:length(headers)
    eval(['kspace.' strrep(headers{i},'.','_') '= list{:,' num2str(i) '};'])
end
kspace.complexdata = [];
kspace.kspace_properties = kspace_properties;
function [nrofheaderlines,kspace_properties] = getNrofheaderlines(listfile)
fid = fopen(listfile);
nrofheaderlines = 0;
currentline = '#';
kspace_properties = struct();
while any(strcmp(currentline(1),{'#','.'}))
    currentline = fgetl(fid);
    nrofheaderlines = nrofheaderlines + 1;
    if strcmp(currentline(1),'.')
        % save property;
        C = textscan(currentline,'. %f %f %f %s : %f %f',1);
        fieldname = strrep(C{4}{1},'-','_');
        kspace_properties.(fieldname) = [C{5} C{6}];
    end
end
nrofheaderlines = nrofheaderlines - 1;
fclose(fid);