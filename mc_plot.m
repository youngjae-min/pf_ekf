
noise = [0 0.125 0.25 0.5 1 1.5 2];
ifagent = 0;

if ifagent
    err_int = mean(err(:,:,:,2:5),4);
    errm = mean(err_int,2);
    errq = quantile(err_int,3,2);
    err_int_h = mean(err_h(:,:,:,2:5),4);
    errm_h = mean(err_int_h,2);
    errq_h = quantile(err_int_h,3,2);
else
    errm = mean(err,2);
    errm_h = mean(err_h,2);
    errq = quantile(err,3,2);
    errq_h = quantile(err_h,3,2);
end



figure
set(gcf,'Units','inches',...
    'Position',[2 2 2.2 1.7],...
    'PaperPositionMode','auto');


for i=3%1:7
    clf
    box on
    grid on
    hold on
    
    p1 = plot(squeeze(errm_h(i,1,:,1)),'color',[0.8500, 0.3250, 0.0980],'linewidth',1.2);
    h = fill([1:100 100:-1:1], [squeeze(errq_h(i,1,:,1))' fliplr(squeeze(errq_h(i,3,:,1))')],'b');
    set(h,'edgecolor','none')
    set(h,'facealpha',.3)
    set(h,'facecolor',[0.9290, 0.6940, 0.1250])
    
    p2 = plot(squeeze(errm(i,1,:,1)),'color',[0 0.6 0.2],'linewidth',1.2);
    o = fill([1:100 100:-1:1], [squeeze(errq(i,1,:,1))' fliplr(squeeze(errq(i,3,:,1))')],'b');
    set(o,'edgecolor','none')
    set(o,'facealpha',.3)
    set(o,'facecolor',[0.4660, 0.6740, 0.1880])
    
    legend([p1 p2],{'PF-only','proposed'})
    
    xlabel('Time [s]')
    ylabel('Distance error [m]')
    axis([0 100 0 20])
    set(gca,...
        'Units','normalized',...
        'YTick',0:5:20,...
        'XTick',0:20:100,...
        'Position',[.11 .16 .855 .82],...
        'FontUnits','points',...
        'FontWeight','bold',...
        'FontSize',7,...
        'FontName','Times')

    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

    print(sprintf('./figures/mc_%03.f.pdf',noise(i)*100),'-dpdf','-painters')
%     print(sprintf('./figures/mca_%03.f.pdf',noise(i)*100),'-dpdf','-painters')
end

% type=1

% for i=1:6
%     subplot(2,6,i)
%     pos = get(gca, 'Position');
%     pos(1) = (i-1)/6+0.01;
%     pos(3) = 1/8;
%     set(gca, 'Position', pos)
%     axis([0 100 0 15])
%     hold on
%     plot(squeeze(errm(i,1,:,type)))
%     plot(squeeze(errm_h(i,1,:,type)))
%     hold off
% end
% 
% for i=1:6
%     subplot(2,6,i+6)
%     pos = get(gca, 'Position');
%     pos(1) = (i-1)/6+0.01;
%     pos(3) = 1/8;
%     set(gca, 'Position', pos)
%     axis([0 100 0 15])
%     hold on
%     plot(squeeze(mean(errm(i,1,:,2:5),4)))
%     plot(squeeze(mean(errm_h(i,1,:,2:5),4)))
%     hold off
% end
