function [kernel,info] = fraunhofer_scattering_kernel(theta,diameter,lambda0,nMedium,variant)
%FRAUNHOFER_SCATTERING_KERNEL Forward diffraction of the projected sphere disk.
% theta rad (0 <= theta < pi/2); lengths m; host index real.
% sphere (default): |x^2*(1+cos(theta))/2 * J1(x*sin(theta))/(x*sin(theta))|^2 / k^2.
% paraxial: sin(theta)->theta and obliquity->1; same units and forward limit.
% This is scalar diffraction, not reflection/refraction or the full Mie field.
    if nargin<4 || isempty(nMedium), nMedium=1; end
    if nargin<5 || isempty(variant), variant='sphere'; end
    variant=validatestring(variant,{'sphere','paraxial'});
    validateattributes(theta,{'numeric'},{'vector','nonempty','real','finite','>=',0,'<',pi/2});
    validateattributes(diameter,{'numeric'},{'vector','nonempty','real','finite','positive'});
    validateattributes(lambda0,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(nMedium,{'numeric'},{'scalar','real','finite','positive'});
    theta=theta(:); diameter=diameter(:).';
    k=2*pi*nMedium/lambda0; x=k*diameter/2;
    if strcmp(variant,'sphere')
        u=sin(theta)*x; obliquity=(1+cos(theta))/2;
    else
        u=theta*x; obliquity=ones(size(theta));
    end
    ratio=zeros(size(u)); small=abs(u)<1e-3;
    % J1(u)/u = 1/2-u^2/16+u^4/384-...; well-defined at theta=0.
    ratio(small)=0.5-u(small).^2/16+u(small).^4/384;
    ratio(~small)=besselj(1,u(~small))./u(~small);
    kernel=(obliquity.*ratio).^2.*(x.^4)/k^2;
    info=struct('model','fraunhofer','variant',variant,'k',k,'x',x, ...
        'nominalDiffractionCrossSection',pi*diameter.^2/4, ...
        'validity','large-particle forward scalar diffraction; not full-sphere scattering');
end
