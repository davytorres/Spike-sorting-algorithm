function highres(filename)
oldscreenunits= get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches','PaperPosition',newpos)
%print('-dtiff',filename,'-r1200');
print('-djpeg',filename,'-r1200');
drawnow
set(gcf,'Units',oldscreenunits,'PaperUnits',oldpaperunits,'PaperPosition',oldpaperpos);
