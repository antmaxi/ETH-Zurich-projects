function printFig(dir, fileName)
    disp(['Printing ', fileName, '...']);
    cmdDir = sprintf('mkdir -p %s', dir);
    system(cmdDir);
    print([dir, fileName],'-dpdf', '-fillpage');
    cmdCrop = sprintf('pdfcrop %s%s %s%s', dir, fileName, dir, fileName);
    [~, ~] = system(cmdCrop);
end