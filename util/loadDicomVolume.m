function [ref, par] = loadDicomVolume(filename, par)

    list = dir(filename);
    
    ref = zeros(par.X_resolution, par.Y_resolution, par.Z_resolution, par.Ne);
    
    for f = 1:numel(list)
        
        info = dicominfo([list(f).folder '/' list(f).name]);
        
        if (info.Width == par.X_resolution)
            
            if ~isfield(par, 'PixelBW')
                par.PixelBW = info.PixelBandwidth;
            end
            
            if ~isfield(par, 'PixelDim')
                par.PixelDim = info.PixelSpacing;
            end
            
            if ~isfield(par, 'TE')
                par.TE = info.EchoTime;
            elseif ~ismember(par.TE, info.EchoTime)
                par.TE = [par.TE info.EchoTime];
            end
            
            if ~isfield(par, 'slicepos')
                par.slicepos = info.ImagePositionPatient';
            elseif ~ismember(par.slicepos, info.ImagePositionPatient', 'rows') 
                par.slicepos = [par.slicepos; info.ImagePositionPatient'];
            end
            
            ref(:, :, ismember(par.slicepos, info.ImagePositionPatient', 'rows'), ...
                par.TE == info.EchoTime) = dicomread([list(f).folder '/' list(f).name]);
            
        end
        
    end

end
