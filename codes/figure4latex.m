function figure4latex(pdf_Name)
DPI = 300; 
%% This will save a "cropped" pdf in the directory
set(gcf,'Units','inches')
h=get(gcf,'Position');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0, 0 ,h(3), h(4)]);
set(gcf, 'PaperSize', [h(3), h(4)])
print('-dpdf',strcat('-r',num2str(DPI)),pdf_Name)