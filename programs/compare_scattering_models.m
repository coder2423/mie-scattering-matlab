function metrics = compare_scattering_models(theta,mieDcs,fraunhoferDcs,forwardMax)
%COMPARE_SCATTERING_MODELS Quantify amplitude, shape and cone-power differences.
% No best-fit scaling; angular curve L2 is integral over theta, not solid angle.
% Cone powers instead use 2*pi*sin(theta)*dtheta (full azimuth average).
    theta=theta(:); mieDcs=mieDcs(:); fraunhoferDcs=fraunhoferDcs(:);
    validateattributes(theta,{'numeric'},{'vector','real','finite','nonnegative'});
    validateattributes(mieDcs,{'numeric'},{'vector','real','finite','nonnegative'});
    validateattributes(fraunhoferDcs,{'numeric'},{'vector','real','finite','nonnegative'});
    if numel(theta)<2 || any(diff(theta)<=0) || theta(1)~=0 || ...
            numel(theta)~=numel(mieDcs) || numel(theta)~=numel(fraunhoferDcs)
        error('Comparison:Grid','Need matching curves on an increasing theta grid starting at zero.');
    end
    validateattributes(forwardMax,{'numeric'},{'scalar','real','finite','positive'});
    if forwardMax>theta(end), error('Comparison:Range','forwardMax exceeds theta grid.'); end
    peakM=max(mieDcs); peakF=max(fraunhoferDcs);
    if peakM<=0
        error('Comparison:ZeroMie','Relative model errors are undefined for a zero Mie curve.');
    end
    difference=fraunhoferDcs-mieDcs;
    forwardTheta=unique([theta(theta<=forwardMax);forwardMax]);
    forwardM=interp1(theta,mieDcs,forwardTheta);
    forwardF=interp1(theta,fraunhoferDcs,forwardTheta);
    shapeF=zeros(size(fraunhoferDcs));
    if peakF>0, shapeF=fraunhoferDcs/peakF; end
    shapeM=mieDcs/peakM;
    mask=mieDcs>=1e-6*peakM;
    solidAngleWeight=2*pi*sin(theta);
    coneM=trapz(theta,mieDcs.*solidAngleWeight);
    coneF=trapz(theta,fraunhoferDcs.*solidAngleWeight);
    forwardConeM=trapz(forwardTheta,forwardM.*(2*pi*sin(forwardTheta)));
    forwardConeF=trapz(forwardTheta,forwardF.*(2*pi*sin(forwardTheta)));
    pointError=max(abs(difference(mask))./mieDcs(mask));
    metrics=struct('angularL2Percent',100*sqrt(trapz(theta,difference.^2)/trapz(theta,mieDcs.^2)), ...
        'forwardAngularL2Percent',100*sqrt(trapz(forwardTheta,(forwardF-forwardM).^2) ...
        /trapz(forwardTheta,forwardM.^2)), ...
        'independentPeakShapeL2Percent',100*sqrt(trapz(theta,(shapeF-shapeM).^2) ...
        /trapz(theta,shapeM.^2)), ...
        'maxPointwisePercentAboveFloor',100*pointError, ...
        'pointwiseFloorFraction',1e-6, ...
        'forwardIntensityDifferencePercent',100*(fraunhoferDcs(1)/mieDcs(1)-1), ...
        'forwardConeDifferencePercent',100*(forwardConeF/forwardConeM-1), ...
        'fullConeDifferencePercent',100*(coneF/coneM-1), ...
        'mieForwardConeCrossSection',forwardConeM,'fraunhoferForwardConeCrossSection',forwardConeF, ...
        'mieFullConeCrossSection',coneM,'fraunhoferFullConeCrossSection',coneF);
end
