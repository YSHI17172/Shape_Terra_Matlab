function GenReport(folderpath)
opts.outputDir = folderpath;
opts.showCode = false;
file = publish('',opts);
web(file)
end