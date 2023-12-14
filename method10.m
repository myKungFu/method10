function [analysis10] = method10(Za,Zb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis10b = method10(Za,Zb);
% analysis10b = method10;
%
% Method10 is an analysis method for determining significance of change.
% This code is an implementation of the methods described in Goodman, Haysely,
% and Jennings (2024), JARO.
%
% Input Arguments:
%  Za = a vector of complex coefficients representing the baseline condition
%  Za = a vector of complex coefficients representing a potential change
%  If the function is called without input arguments, simulated data will
%  be generated and used.
%
% Output Argument:
%  analysis10 = a structure containing the results of the analysis
%
% Author: Shawn Goodman, PhD
%         Dept. of Communication Sciences and Disorders
%         The University of Iowa
% Date: Decembrer 13, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doPlot = 0; % switch to turn on and off plotting

    if nargin == 0
        % create simulated data set
        changeType = 'both'; % 'mag' or 'phi' or 'both'
                            % direction will be picked randomly
                            % (increase/decrease or lead/lag)
        tqc = .25; % total quanity of change must be between 0 and 1, inclusive (0 and 6 dB)
        snr = 30; % signal to noise ratio of simulated dat
        N = 250; % number of sweeps; the length of Za and Zb. 
        [Za,Zb,SNR,Usc] = stimQAD(N,snr,tqc,changeType); % simulate the data
    end

    Za = Za(:); % force data to be column vectors.
    Zb = Zb(:);

    % rotate and scale the matrices
    theta = angle(sum(Za));
    k = abs(mean(Za));
    Za = Za * (1/k)*exp(-1i*theta);
    Zb = Zb * (1/k)*exp(-1i*theta);
    Delta = Zb-Za+1; % complex ratio
    
    % reject extreme values
    gm = geoMed([real(Delta),imag(Delta)]);
    gm = gm(1) +1i*gm(2); % geometric mean in complex form
    ed = abs(Delta - gm); % Euclidian distance from the geometric median
    ed = ed(:)';
    indx = AR(ed,'mild',0);
    nRejectsD = length(indx);
    pcExtreme = nRejectsD / length(Za); % percentage of extreme values
    % remove the rejections calculated by the difference
    Za = AR_engine(Za.',indx);
    Za = Za(:);
    Zb = AR_engine(Zb.',indx);
    Zb = Zb(:);

    % recalculate the complex ratio
    theta = angle(sum(Za));
    k = abs(mean(Za));
    Za = Za * (1/k)*exp(-1i*theta);
    Zb = Zb * (1/k)*exp(-1i*theta);
    Delta = Zb-Za+1; % complex ratio

    [~,~,snrA] = getSEM(Za); % calcualte SNR of the baseline
    tdMin = 3.107*exp(-0.1136*snrA); % minimum detectable turned delta, given snr
    magIncMin = 20*log10(1+tdMin);
    if tdMin > 1
        magDecMin = -12;
    else
        magDecMin = 20*log10(1-tdMin);
    end
    phiMin = abs(2*asin(tdMin/2));


    [lambda,upsilon,g,R,delta] = toleranceRegion(Delta);  % Compute statistical tolerance region
    [sig,sigMag,sigPhi,sigBoth1,sigBoth2] = determineSignificance(lambda,upsilon,g,R);

    td = abs(delta-1); % turned delta (total change)
    deltaMag = abs(delta); % change in magnitude
    deltaPhi = angle(delta); % change in phase
    pcMag = abs(deltaMag-1) / td; % percentage of change due to magnitude
    pcPhi = 1 - pcMag; % percentage of change due to phase

    % save output in a structure
    analysis10.delta = delta;
    analysis10.td = td;
    analysis10.deltaMag = deltaMag;
    analysis10.deltaPhi = deltaPhi;
    analysis10.pcMag = pcMag;
    analysis10.pcPhi = pcPhi;
    analysis10.sig = sig;
    analysis10.sigMag = sigMag;
    analysis10.sigPhi = sigPhi;
    analysis10.sigBoth1 = sigBoth1;
    analysis10.sigBoth2 = sigBoth2;
    analysis10.snrA = snrA;
    analysis10.pcExtreme = pcExtreme;
    analysis10.tdMin = tdMin;
    analysis10.magIncMin = magIncMin;
    analysis10.magDecMin = magDecMin;
    analysis10.phiMin = phiMin;

    if doPlot == 1
        Theta = linspace(-pi,pi,4000)';
        U =  exp(1i*Theta);
        figure
        plot(U,'k')
        hold on
        plot(Delta,'.k')
        plot(R,'r')
        plot(delta,'r*')
        plot(Usc,'g')
        grid on
        line([-1 1],[0 0],'Color',[0 0 0],'LineStyle','-','LineWidth',0.5)
        line([0 0],[-1 1],'Color',[0 0 0],'LineStyle','-','LineWidth',0.5)
        title(['Sig=',num2str(sig),'  SigMag=',num2str(sigMag),'   SigPhi=',num2str(sigPhi),'  SigBoth1=',num2str(sigBoth1),'   SigBoth2=',num2str(sigBoth2)])    
    end
end
% INTERNAL FUNCTIONS ------------------------------------------------------
function [lambda,upsilon,g,R,delta] = toleranceRegion(Delta)
    % Calculate the ellipse R describing the standard error of the mean 
    % around the data in vector D.
    Delta = Delta(:); % force to a column vector
    M = length(Delta);
    delta = mean(Delta,1); % mean of the bootstrap estimates

    Delta = Delta - delta; % remove the mean to center around zero
    Delta = [real(Delta),imag(Delta)]; % convert to two-column matrix
    S = (1/(M-1))*(Delta'*Delta); % variance, covariance matrix
    [upsilon,lambda] = eig(S); 
    % upsilon is a 2x2 matrix with columns of eigen vectors
    % lambda is a diagonal matrix of eigenvalues

    m = 2000; % number of samples in the ellipse
    theta = linspace(0,2*pi,m); % Linearly-spaced vector 
    C = [cos(theta);sin(theta)]; % Create unit circle
    g = 0.385*sqrt(M); % scaling factor.

    R = (upsilon*sqrt(lambda)*C)' / g; % elliptical region
    R = R(:,1) + 1i*R(:,2); % convert back to complex vector
    R = R + delta; % put back to the original location around delta

    lambda = diag(lambda); % just save the diagonal values (off-diagonal values are zero)
end
function [SIG,MAG,PHI,BOTH1,BOTH2] = determineSignificance(lambda,upsilon,g,R)
    % DETERMINE STATISTICAL SIGNIFICANCE
    % sig = overall significance
    % MAG = 1 = mag is significant, either alone or with phase
    % PHI = 1 = phase is significant, either alone or with phase
    % BOTH1 = 1 = mag alone and phase alone are both significant
    % BOTH2 = 1 = neither mag alone nor phase alone are significant, but
    %             together they are significant
    %
    % ALGORITHM:
    % Determine whether R intersects
    %   1) The point 1+0i (if so, not significant)
    %   2) the unite circle anywhere (if so, magnitude not significant) 
    %   3) the x-axis (if so, phase not significant)

    % Check whether the imaginary part of R contains zero
    % If it does, phase alone is not significant
    im0 = (min(imag(R))<=0 & max(imag(R))>=0); 

    % Make R a standard ellipse centered at origin and rotated so 
    % major axis lies on the x-axis
    smallUpsi = min(upsilon(1:2)); 
    bigUpsi = max(upsilon(1:2));
    theta = atan2(bigUpsi,smallUpsi); % amount by which the ellipse is rotated
    
    % create a line segment for U (unit circle);
    Theta = linspace(-pi,pi,4000)';
    U =  exp(1i*Theta);

    % shift the origin to the center of the ellipse
    shift = mean(R);
    R = R - shift;
    U = U - shift;
    O = complex(1+1i*0) - shift; % the "no change" point 1 + 0i
    % rotate to make the ellipse standard
    R = R * exp(1i*theta);
    U = U * exp(1i*theta);
    O = O * exp(1i*theta);
    % find a and b values (lengths of major and minor axes) of R
    ab = sqrt(lambda)/g;
    a = max(ab); % length of major axis
    b = min(ab); % length of minor axis
    
    % For each point, solve p = (x^2/a^2) + (y^2/b^2);
    % If p is < 1 the point is inside the ellipse
    % If p = 1, the point is on the ellipse
    % If p > 1, the point is outside the ellipse
    
    % Overall significance: Does the ellipse contain the point 1 + 0i? (O)
    p0 = (real(O).^2./a^2) + (imag(O).^2./b^2);

    % Is any of the unit circle inside R? If so, Mag alone is not significant
    p = real(U).^2./a^2 + imag(U).^2./b^2;
    inIndx = find(p<1); 
    
    % Calculate significance:
    %   isempty(inIndx) means mag is significant
    %   im0 = 0 means phase is significant
    %   p0 <1 means origin inside (sig = 0)
    
    if isempty(inIndx)
        MAG = 1;
    else
        MAG = 0;
    end
    if im0 == 0
        PHI = 1;
    else
        PHI = 0;
    end
    if p0<1
        SIG = 0;
    else
        SIG = 1;
    end
    if MAG==0 && PHI==0 && SIG==1
        BOTH2 = 1;
    else
        BOTH2 = 0;
    end
    if MAG==1 && PHI==1
        BOTH1 = 1;
        % by leaving this part off, easily can ask if mag is significant
        % (either alone or with phase) and if phi is significant (either
        % alone or with phase), which is how we will use in the paper.
        %MAG = 0;
        %PHI = 0;
    else
        BOTH1 = 0;
    end

end
function [Za,Zb,SNR,U] = stimQAD(N,snr,tqc,changeType)
    % make stimulus (quick and dirty method, i.e., not using sinusoids)
    noiseLvl_dB = -40; % noise floor level
    ref = 0.000001; % reference in microvolts
    noiseLvl = 10^(noiseLvl_dB/20)*ref; % noise level
    offset = 21;
    sigLvl_dB = noiseLvl_dB + snr - offset; % signal level
    sigLvl = 10^(sigLvl_dB/20)*ref; % noise level

    % create without CAS (A)
    mag = sigLvl;
    phi = 0;
    x = mag*cos(phi);
    y = mag*sin(phi);
    z = complex(x + 1i*y);
    noise_x = randn(N,1)*noiseLvl;
    noise_y = randn(N,1)*noiseLvl;
    noise_z = complex(noise_x + 1i*noise_y);
    Za = z + noise_z;

    K = length(Za);
    Xk = Za;
    Xbar = (1/K) * sum(Xk); % signal is the coherent mean
    Xbar2 = abs(Xbar) .^2; % signal energy
    XBAR = Xbar; % replicate to matrix (to avoid running a loop)
    S2 = (1/(K-1)) * sum((Xk - XBAR) .* conj(Xk - XBAR)); % variance
    Se2 = (1/K) * S2; % energy of the standard error
    SNR = 10*log10(Xbar2/Se2); % signal to noise ratio of without CAS

    % create with CAS (B)
    % here, create an equivalent change (must be 6 dB or less) of either
    % magnitude, phase, or some random combination of the two
    %
    % tqc = abs(delta-1)
    % mag increase delta = 1+tqc
    % mag decrease = 1-tqc
    % phi = abs(2*asin((1-abs(1+tqc))/2))
    %tqc = [0, .25, .333, .5, .666, .75, 1]'; % tqc is percentage change, so max is 100% increase = 6 dB change
    %deltaInc = 1 + tqc;
    %deltaDec = 1 - tqc;
    %phi = abs(2*asin((1-abs(1+tqc))/2)); % max phase change is pi/3

    Theta = linspace(-pi,pi,4000)';
    U =  exp(1i*Theta);
    U = U * tqc; % scale unit circle by magnitude
    U = U + 1; % shift to the point 1+0i


    % randomly choose increase or decrease (phase lead or lag)
    coin = randn(1);
    if coin <= 0
        direction = 'dec';
    else
        direction = 'inc';
    end

    if strcmp(changeType,'mag')
        if strcmp(direction,'inc')
            dbchange = 20*log10(1+tqc);
        elseif strcmp(direction,'dec')
            dbchange = 20*log10(1-tqc);
        end
        sigLvl_dB = noiseLvl_dB + snr - offset + dbchange; % signal level
        sigLvlB = 10^(sigLvl_dB/20)*ref; % noise level
        mag = sigLvlB;
        phi = 0;
    elseif strcmp(changeType,'phi')
        dbchange = 20*log10(1+tqc);
        mag = sigLvl;
        x = abs(dbchange);
        x = 10^(x/20); % dB change expressed as a ratio
        phi = abs(2*asin((1-abs(x))/2)); % equivalent phase change
        if strcmp(direction,'dec')
            phi = -phi;
        end
    elseif strcmp(changeType,'both')
        % in this case, randomly choose a combintation of mag and phase 
        indx = floor(rand(1)*(length(U)));
        draw = U(indx);
        mag = abs(draw);
        phi = angle(draw);
        mag = mag * sigLvl;
    end
    x = mag*cos(phi);
    y = mag*sin(phi);
    z = complex(x + 1i*y);
    noise_x = randn(N,1)*noiseLvl;
    noise_y = randn(N,1)*noiseLvl;
    noise_z = complex(noise_x + 1i*noise_y);
    Zb = ones(size(noise_z))*z + noise_z;

    %plot(Za,'.r')
    %hold on
    %plot(Zb,'.b')
end
function [bar,sem,snr] = getSEM(FD)
    bar = abs(mean(FD));
    sem = std(FD)/sqrt(length(FD));
    snr = 20*log10(bar/sem);
    % alterative calculation:
    %K = length(FD);
    %bar = mean(FD); % coherent complex mean
    %S2 = (1/(K-1)) * sum((FD - bar) .* conj(FD - bar)); % variance
    %sem = sqrt((1/K) * S2); % standard error of mean
    %snr = 20*log10(abs(bar)/sem);
end
function [GM] = geoMed(X)
    % compute the geometric median by iteratively reducing the error
    W=ones(size(X,1),1); 
    W=W/sum(W);
    GM = W'*X;
    GM = GM(:)';
    epsilon = inf; % error term
    
    regulate = 0.0001;
    maxIterations = 100;
    minError = 1E-7;
    iteration = 0;
    while (iteration<=maxIterations) && (epsilon>minError)
        iteration = iteration + 1;
        w = sqrt(sum(X.^2,2))./W;
        w = 1./(w+0.00001); 
        GMi = sum(X.*w,1)/sum(w);
        epsilon = sum((GM - GMi).^2);
        GM=GMi;
        regulate = max(regulate/10,eps); 
    end
end
function [indx,nRejects] = AR(y,tolerance,doPlot)
    y = y(:); % ensure y is a column vector
    warning('OFF')
    if nargin <= 2
        doPlot = 0; % set = 1 to look at outlier plot
    end
    if nargin == 1,
        tolerance = 'moderate';
    end
    switch tolerance
        case 'mild'
            multiplier = 1.5;
        case 'moderate'
            multiplier = 2.25;
        case 'severe'
            multiplier = 3;
        otherwise
            disp('Error: requested tolerance is not correctly specified.')
            return
    end
    
    xx = sort(y(:));
    [dummy,indx] = max(xx);
    xx = xx(1:indx); % get rid of any nan
    N = length(xx); % number of observations
    q = 100 *(0.5:N-0.5)./N;
    xx = [min(xx); xx(:); max(xx)];
    q = [0 q 100];
    F1 = interp1(q,xx,25); % the first fourth, approx 25th precentile
    F2 = interp1(q,xx,50); % the first half, approx 50th precentile
    F3 = interp1(q,xx,75); % the third fourth, approx 75th percentile
    IQR = F3 - F1; % the interquartile range
    
    ArtifactVector = y >= (F1-multiplier*IQR) & y <= (F3 + multiplier*IQR);
    [indx,val] = find(~ArtifactVector); % index the artifacts which should be rejected
    
    nRejects = length(indx); % number of rejected buffers
    percentRejected = (length(indx) / N) * 100;
    % disp([num2str(percentRejected),' percent of the buffers rejected as artifacts.'])
    % commented out previous line on 6/3/2010--ras and ssg
    
    if doPlot == 1,
        figure(63)
        plot(q,xx,'*b-')
        hold on
        L1 = line([q(1) q(end)],[F1 F1]);
        L2 = line([q(1) q(end)],[F2 F2]);
        L3 = line([q(1) q(end)],[F3 F3]);
        L4 = line([q(1) q(end)],[F1-multiplier*IQR F1-multiplier*IQR]);
        L5 = line([q(1) q(end)],[F3 + multiplier*IQR F3 + multiplier*IQR]);
        set(L1,'Color',[0 1 0])
        set(L2,'Color',[0 0 0],'LineStyle',':')
        set(L3,'Color',[0 1 0])
        set(L4,'Color',[1 0 0])
        set(L5,'Color',[1 0 0])
        hold off
        %pause
    end
    warning('ON')
end
function [Y] = AR_engine(Y,indx)
    % This function removes buffers with artifacts from a data set.
    Y(:,indx) = [];
end