function [Gene, Gene2] = shownameskat(FILE,GENES)

% NAME = shownames(FILE,I,GENES)
% This function determines the names of the genes
% FILE is the filename of the text file that contains all the genenames
% I is the index matrix containing the gene numbers (row numbers) for the gene names that have to be determined
% GENES is the total number of genes in the text file

V = fopen(['d:\matlab\testen\yeast130601\',FILE],'r');
i=0;
for row = 1:GENES
   line = fgetl(V);
   Gene(row,1:size(line,2))=line;
   Gene2{row,1}=line;


end

fclose('all');


