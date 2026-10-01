function c=validate_scattering_config(c)
%VALIDATE_SCATTERING_CONFIG Check settings before expensive kernel construction.
    if ~isstruct(c) || ~isscalar(c), error('Simulation:Config','Configuration must be a scalar struct.'); end
    validateattributes(c.wavelengths_nm,{'numeric'},{'vector','nonempty','real','finite','positive'});
    validateattributes(c.diameterRange_um,{'numeric'},{'vector','numel',2,'real','finite','positive'});
    validateattributes(c.angleRange_deg,{'numeric'},{'vector','numel',2,'real','finite','nonnegative','<',90});
    validateattributes(c.diameterPoints,{'numeric'},{'scalar','integer','>=',2});
    validateattributes(c.anglePoints,{'numeric'},{'scalar','integer','>=',3});
    validateattributes(c.forwardMax_deg,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(c.nMedium,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(c.nParticle,{'numeric'},{'scalar','finite'});
    if real(c.nParticle)<=0 || imag(c.nParticle)>0
        error('Mie:IndexConvention','Use a passive absolute particle index n-i*kappa with n>0.');
    end
    if c.nParticle==c.nMedium
        error('Comparison:ZeroMie','Index-matched particles have zero Mie scattering; relative metrics are undefined. Use the single-sphere or distribution API instead.');
    end
    validateattributes(c.incidentIrradiance,{'numeric'},{'scalar','real','finite','nonnegative'});
    validateattributes(c.distance_m,{'numeric'},{'scalar','real','finite','positive'});
    if ~isempty(c.totalNumber), validateattributes(c.totalNumber,{'numeric'},{'scalar','real','finite','nonnegative'}); end
    c.polarization=validatestring(c.polarization,{'unpolarized','parallel','perpendicular','linear'});
    validateattributes(c.azimuth_rad,{'numeric'},{'real','finite','nonempty'});
    if ~isscalar(c.azimuth_rad) && (~isvector(c.azimuth_rad) || numel(c.azimuth_rad)~=c.anglePoints)
        error('Mie:Azimuth','Azimuth must be scalar or match anglePoints.');
    end
    c.fraunhoferVariant=validatestring(c.fraunhoferVariant,{'sphere','paraxial'});
    if numel(unique(c.wavelengths_nm))~=numel(c.wavelengths_nm)
        error('Simulation:Wavelengths','Wavelengths must be distinct to avoid duplicate output files.');
    end
    if ~(ischar(c.outputFolder) && isrow(c.outputFolder)) && ~(isstring(c.outputFolder) && isscalar(c.outputFolder) && ~ismissing(c.outputFolder))
        error('Simulation:OutputFolder','outputFolder must be a nonempty text scalar.');
    end
    c.outputFolder=char(c.outputFolder);
    if isempty(strtrim(c.outputFolder)), error('Simulation:OutputFolder','outputFolder must not be empty.'); end
    validateattributes(c.makePlots,{'logical'},{'scalar'}); validateattributes(c.saveData,{'logical'},{'scalar'});
    validateattributes(c.saveFigures,{'logical'},{'scalar'});
    if c.saveFigures && ~c.makePlots
        error('Simulation:FigureSaving','Set makePlots=true to save figures; use figureVisible=off for background export.');
    end
    if c.diameterRange_um(2)<=c.diameterRange_um(1) || c.angleRange_deg(1)~=0 || c.angleRange_deg(2)<c.forwardMax_deg
        error('Simulation:Grid','Use increasing size range and an angle range starting at 0 and covering forwardMax.');
    end
    c.distributionBasis=validatestring(c.distributionBasis,{'volume','number'});
    xMin=pi*c.nMedium*c.diameterRange_um(1)*1e3/max(c.wavelengths_nm);
    xMax=pi*c.nMedium*c.diameterRange_um(2)*1e3/min(c.wavelengths_nm);
    if xMin<1e-6 || xMax>1e4, error('Mie:SizeRange','All diameter/wavelength combinations must satisfy 1e-6 <= x <= 1e4.'); end
    c.plotQuantity=validatestring(c.plotQuantity,{'own_peak','common_peak','dcs','irradiance'});
    c.intensityYScale=validatestring(c.intensityYScale,{'linear','log'});
    c.figureVisible=validatestring(c.figureVisible,{'on','off'});
    if ~(iscellstr(c.distributionTypes) || isstring(c.distributionTypes) || ischar(c.distributionTypes))
        error('Simulation:Types','distributionTypes must contain distribution names.');
    end
    types=cellstr(lower(string(c.distributionTypes)));
    if isempty(types) || numel(unique(types))~=numel(types), error('Simulation:Types','Choose distinct distribution types.'); end
    if ~all(ismember(types,{'rr','normal','lognormal','bimodal'})), error('Simulation:Types','Supported types: rr, normal, lognormal, bimodal.'); end
    c.distributionTypes=types(:).';
    if ~isstruct(c.parameters) || ~isscalar(c.parameters) || ~all(isfield(c.parameters,types))
        error('Simulation:Distribution','Supply parameter structs for every selected distribution.');
    end
    if strcmp(c.plotQuantity,'irradiance') && isempty(c.totalNumber)
        error('Simulation:Irradiance','Set totalNumber before requesting an irradiance plot.');
    end
    validateattributes(c.showPeakDiagnostics,{'logical'},{'scalar'});
    validateattributes(c.referenceDiameter_um,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(c.diagnosticStd_um,{'numeric'},{'vector','numel',2,'real','finite','positive'});
    validateattributes(c.sideLobeYMax,{'numeric'},{'scalar','real','finite','positive','<',1});
    if c.showPeakDiagnostics && (c.referenceDiameter_um<c.diameterRange_um(1) || ...
            c.referenceDiameter_um>c.diameterRange_um(2))
        error('Simulation:ReferenceDiameter','Diagnostic reference diameter must lie within the size range.');
    end
    if c.showPeakDiagnostics && diff(c.diameterRange_um)/(c.diameterPoints-1)>min(c.diagnosticStd_um)/3
        error('Simulation:DiagnosticGrid','Increase diameterPoints or diagnosticStd_um: need >=3 intervals per diagnostic standard deviation.');
    end
end
