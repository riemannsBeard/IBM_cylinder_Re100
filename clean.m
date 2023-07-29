function [] = clean()

clc
clear global
close all


if exist('./forces', 'dir')
    rmdir('./forces', 's')
end

if exist('./stored', 'dir')
    rmdir('./stored', 's')
end

if ~isempty(dir('*.mat'))
    delete *.mat
end

end

