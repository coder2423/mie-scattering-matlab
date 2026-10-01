function plot_scattering_comparison(s,w,targetAxes,view)
%PLOT_SCATTERING_COMPARISON Plot existing results in an array of axes.
% targetAxes: nDistributions-by-3; PSD, forward zoom, and configured angle range.
% view controls display only; no simulation data or kernels are changed.
    c=s.config;
    if nargin<4, view=c; end
    nTypes=numel(s.distributions);
    if ~isequal(size(targetAxes),[nTypes,3]) || ~all(isgraphics(targetAxes(:)))
        error('Simulation:PlotAxes','Provide one row of three axes for each distribution.');
    end
    for j=1:nTypes
        ax=targetAxes(j,1); axZoom=targetAxes(j,2); axFull=targetAxes(j,3);
        for a=targetAxes(j,:), cla(a,'reset'); end
        plot(ax,s.diameter*1e6,s.distributions(j).density/1e6, ...
            'Color',[0.12,0.35,0.65],'LineWidth',1.6);
        xlabel(ax,'Diameter (\mum)'); ylabel(ax,' density (\mum^{-1})');
        if strcmp(c.distributionBasis,'volume'), basisLabel='Volume'; else, basisLabel='Number'; end
        title(ax,[s.distributions(j).name,': ',basisLabel,' PSD']);
        xlim(ax,c.diameterRange_um); grid(ax,'on');
        [yM,yF,yLabel]=scattering_plot_values(s.mieResults{j,w},s.fraunhoferResults{j,w},view.plotQuantity);
        hold(axZoom,'on');
        plot(axZoom,s.thetaDeg,yM,'LineWidth',1.5,'Color',[0.12,0.35,0.65]);
        plot(axZoom,s.thetaDeg,yF,'--','LineWidth',1.5,'Color',[0.85,0.33,0.10]);
        xlim(axZoom,[0,c.forwardMax_deg]); grid(axZoom,'on');
        xlabel(axZoom,'Scattering angle (deg)'); ylabel(axZoom,yLabel);
        title(axZoom,sprintf('Forward 0-%g deg',c.forwardMax_deg));
        legend(axZoom,{'Mie','Fraunhofer'},'Location','northeast'); hold(axZoom,'off');
        hold(axFull,'on');
        plot(axFull,s.thetaDeg,yM,'LineWidth',1.5,'Color',[0.12,0.35,0.65]);
        plot(axFull,s.thetaDeg,yF,'--','LineWidth',1.5,'Color',[0.85,0.33,0.10]);
        set(axFull,'YScale',view.intensityYScale); xlim(axFull,c.angleRange_deg); grid(axFull,'on');
        xlabel(axFull,'Scattering angle (deg)'); ylabel(axFull,yLabel);
        rowIndex=(w-1)*nTypes+j;
        if strcmp(view.plotQuantity,'own_peak')
            title(axFull,sprintf('Normalized shape difference %.2f%%',s.metrics.independentPeakShapeL2Percent(rowIndex)));
        else
            title(axFull,sprintf('Physical curve L2 difference %.2f%%',s.metrics.angularL2Percent(rowIndex)));
        end
        hold(axFull,'off');
        set(targetAxes(j,:),'FontName','Helvetica','FontSize',10,'XScale','linear');
        for a=targetAxes(j,:), a.Toolbar.Visible='off'; end
    end
end
