function FoInfo = CheckFolder(filename)

%CommonFileName = 'x3175.txt';  %or whatever the common name is
FoInfo = dir('..\Output');
FoInfo(~[FoInfo.isdir]) = []; % Remove non-folder entries in folder info
RemIdx(1) = find(strcmp({FoInfo.name},{'.'})==1);
RemIdx(2) = find(strcmp({FoInfo.name},{'..'})==1);
FoInfo(sort(RemIdx)) = []; % Remove . and .. entries from folder info
for FoIdx = length(FoInfo):-1:1
  %specificname = fullfile(FoInfo(FoIdx}.name, filename);
%   if exists(specificname, 'file)
%     %at this point, insert your code to examine specificfile
%   end
    if strfind(FoInfo(FoIdx).name,filename)
        disp(['Match found with ' filename ': ' FoInfo(FoIdx).name]);
        dir(['..\Output\' FoInfo(FoIdx).name])
    end
end

end