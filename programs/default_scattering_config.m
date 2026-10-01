function config = default_scattering_config()
%DEFAULT_SCATTERING_CONFIG Editable defaults for the standalone simulation.
% Range units are indicated by field names. All PSD length parameters are m.
    folder=fileparts(fileparts(mfilename('fullpath')));

    % -------------------- Optical conditions --------------------
    config.wavelengths_nm=[450,632];    % Vacuum wavelengths, nm; one run per wavelength.
    config.nParticle=1.591;             % Absolute index; passive absorption: 1.591-0.01i.
    config.nMedium=1.0;                 % Real index of a nonabsorbing host; air by default.

    % -------------------- Diameter and angle grids --------------------
    config.diameterRange_um=[10,200];   % Diameter interval, um; PSDs are truncated here.
    config.diameterPoints=401;          % Quadrature nodes; refine for narrow distributions.
    config.angleRange_deg=[0,10];       % Forward interval, deg; starts at 0, ends <90.
    config.anglePoints=2001;            % Linear sampling; default spacing is 0.005 deg.
    config.forwardMax_deg=1;            % Forward plot and cone-integral upper limit, deg.

    % -------------------- Particle size distributions --------------------
    config.distributionTypes={'rr','normal','lognormal','bimodal'};
    config.distributionBasis='volume';  % 'volume' or 'number'; weighting is explicit.
    config.parameters.rr=struct('Mode',100e-6,'Shape',3); % Mode requires Shape>1.
    config.parameters.normal=struct('Mean',100e-6,'Std',20e-6);
    config.parameters.lognormal=struct('Mode',100e-6,'LogSigma',0.30); % SD of ln(d).
    config.parameters.bimodal=struct('Modes',[80,130]*1e-6, ...
        'LogSigmas',[0.15,0.10],'VolumeFractions',[0.5,0.5]);
    % Fractions apply within the truncated interval. For number-basis inputs
    % the legacy VolumeFractions field denotes number shares, not volume.

    % -------------------- Polarization and absolute scaling --------------------
    config.polarization='unpolarized';  % Also: parallel, perpendicular, linear.
    config.azimuth_rad=0;               % For linear: angle of incident E to scattering plane.
    config.incidentIrradiance=1;         % W/m^2 at the sample; not laser power in W.
    config.distance_m=1;                % Far-field distance, m; not focal length.
    config.totalNumber=[];              % Total illuminated count; empty disables irradiance.
    config.fraunhoferVariant='sphere';   % sin(theta) and obliquity; 'paraxial' uses theta.
    % Fraunhofer is scalar; polarized Mie does not add polarization to diffraction.

    % -------------------- Plot quantities and peak diagnostics --------------------
    config.plotQuantity='own_peak';     % own_peak, common_peak, dcs, or irradiance.
    config.intensityYScale='log';        % Tail plot: log or linear; x axes stay linear.
    config.showPeakDiagnostics=true;
    config.referenceDiameter_um=100;    % Single-size diameter and normal-PSD center, um.
    config.diagnosticStd_um=[2,20];      % Narrow and broad diagnostic widths, um.
    config.sideLobeYMax=0.04;            % Zoom of the SAME data; clips main-peak top.

    % -------------------- Preview and independent save switches --------------------
    config.makePlots=true;              % false: compute without opening figures.
    config.figureVisible='on';          % 'off' enables invisible batch plots.
    config.saveData=false;              % true: CSV and MAT; default is preview only.
    config.saveFigures=false;           % true: PNG and FIG; independent of saveData.
    config.outputFolder=fullfile(folder,'results','model_comparison');
    % Explicit saving overwrites files with matching names in outputFolder.
end
