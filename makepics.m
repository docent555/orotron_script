function makepics(RES)

% v = VideoWriter('movie_test.avi');
% open(v);

XMAX = max(RES.ZAxis);
YMAXB = max(max(abs(RES.OUTB))*1.1);
YMINB = min(min(abs(RES.OUTB))*1.1);
YMAXJ = max(max(abs(RES.OUTJ))*1.1);
YMINJ = min(min(abs(RES.OUTJ))*1.1);
 
FigHandle = figure();

FolderB = 'pics/B/';
FolderJ = 'pics/J/';

% hash = datestr(now,30);
% FolderName=sprintf('%s%s', Folder, hash);
FolderNameB=sprintf('%s', FolderB);
mkdir(FolderNameB);
FolderNameJ=sprintf('%s', FolderJ);
mkdir(FolderNameJ);

for k=1:length(RES.TAxis(:,1))        
       
    fileB = sprintf('%s/%s', FolderNameB, sprintf('%08.4f.bmp', RES.TAxis(k)));
    
    plot(RES.ZAxis(:,1), abs(RES.OUTB(:,k))); 
    title(sprintf('Time: %.4f [ns]' , RES.TAxis(k,1)));       
    
    xlim([0 XMAX]);
    ylim([YMINB YMAXB]);
    xlabel('Z axis [mm]', 'fontsize', 12);
    
%     frame = getframe(gcf);
%     writeVideo(v,frame);
    saveas(FigHandle, fileB);
    
    fileJ = sprintf('%s/%s', FolderNameJ, sprintf('%08.4f.bmp', RES.TAxis(k)));
    
    plot(RES.ZAxis(:,1), abs(RES.OUTJ(:,k))); 
    title(sprintf('Time: %.4f [ns]' , RES.TAxis(k,1)));       
    
    xlim([0 XMAX]);
    ylim([YMINJ YMAXJ]);
    xlabel('Z axis [mm]', 'fontsize', 12);
    
%     frame = getframe(gcf);
%     writeVideo(v,frame);
    saveas(FigHandle, fileJ);
    
    disp(k)
    
end

% close(v);
end

