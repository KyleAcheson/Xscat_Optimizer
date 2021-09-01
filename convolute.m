function [conv_signal, tconv, width] = convolute(signal, time, fwhm, debug, qAng)

    % TO-DO: Merged convolution of individual pump and probe pulses
    % KA - 24/5/2021 - revised to use circulant matrices - v fast
    % Convolute a signal with a single gaussian pulse with a set fwhm
    % (normalised to unity)
    % signal is padded by 3*fwhm before t=0 and after t=end
    
    % INPUTS:
    % signal - original input signal - time automatically set to column
    % time - original time vector for the unconvoluted signal
    % fwhm - fwhm of pulse - wil combine pump and probe fwhm together
    % debug - 0 = no debugging, 1 = time convolution and inspect kernel
    
    % OUTPUTS:
    % conv_signal - convoluted signal with same dimensions as original
    % signal. The padding is cut off before being returned.

    if debug == 1 tstart = tic; end
    
    fwhm = sum(fwhm); % sum of two pump and probe pulse fwhms
    if fwhm == 0; warning("FWHM set to zero - no convolution done."); end

    [nr, nc] = size(signal); % size input signal
    nt = length(time);
    difftime = diff(time);
    dt = sum(difftime)/(nt-1); % time step
    width = fwhm/(2*sqrt(2*log(2)));
    duration = ceil(3*width); % length by which to pad 
    

    if dt ~= difftime(1) warning("Non-linear spacing in time."); end
    if nr ~= nt && nc ~= nt error("Inconconsistent dimensions."); end
    if nc ~= nt signal = signal.'; [nr, nc] = swap(nc, nr); end % ensure time always 2nd dim
   

    tmin = min(time); tmax = max(time);
    multiple = 50;
    tconv_min = tmin - duration ; tconv_max = tmax + duration; % pad by 3*sigma
    tconv_min = roundn(tconv_min, multiple); tconv_max = roundn(tconv_max, multiple);
    
    tconv = tconv_min:dt:tconv_max; % extended convolution time
    last_tconv = find(tconv == time(end)); % last point where full width of gaussian present
    ntc = length(tconv); 
    padend = round(duration/dt); % new time zero index
    padstart = round(padend + nt-1); % new 'actual' signal end index
    
    
    
    tmat = tmin-(duration*2):dt:tmax+(duration*2);
    Kvec = gaussian(fwhm, tmat, tconv(1));
    Kernel =  toeplitz([Kvec(1) fliplr(Kvec(2:end))], Kvec);
  
    inds = ismember(tmat, tconv);
    Kernel = Kernel(1:ntc, inds);
  
    %Kvec = gaussian(fwhm, tconv, tconv(1)); % init centre of gaussian on t=0
    %Kernel = toeplitz([Kvec(1) fliplr(Kvec(2:end))], Kvec); % circulant matrix - like an integral kernel
    
    pad_signal = zeros(nr, ntc);

    for i=1:padend
    pad_signal(1:nr, i) = signal(1:nr, 1); % t < 0 - pad by setting to t=0 value (extended by length 3*sigma)
    end
    
    %keyboard
    pad_signal(1:nr, padend:padstart) = signal; % original signal in middle
    
    for i=padstart:ntc
    pad_signal(1:nr, i) = signal(1:nr, end); % t > original end time - pad by setting to t=end value 
    end
    
    if debug == 1
        [QQ, TT] = meshgrid(qAng, tconv);
        figure
        mesh(QQ, TT, pad_signal.')
        title('Padded signal pre-convolution')
        keyboard
    end

    conv_signal = zeros(nr, ntc);

    conv_signal = pad_signal * Kernel'; % perform convolution as matmul w gaussian integral kernel
    
  
    conv_signal = conv_signal(:, 1:last_tconv); % return convoluted signal with original dims
    tconv = tconv(1:last_tconv);
    
    if all(isnan(conv_signal(1,:))) == 1; conv_signal(isnan(conv_signal)) = 0; warning("NaNs at time zero set to 0."); end

    if debug == 1
        
        [TT, TT] = meshgrid(tconv, tconv);
        [QQ, TQ] = meshgrid(qAng, tconv);
        tind = find(tconv == 0);
        %keyboard
        figure
        subplot(1,2,1)
        plot(tconv, Kernel(tind, 1:last_tconv))
        xticks(-500:100:1500)
        title('Inspection of Kernel at original t=0')
        subplot(1,2,2)
        imagesc(Kernel) % 2D view of Kernel
        title(['2D view of Kernel. Time zero index =', num2str(tind)])
        disp(['Convolution integrates to: ', num2str(sum(Kernel(1,:)))])
        telapsed = toc(tstart);
        disp(['Time elapsed for CONV   (s):' num2str(telapsed)])
        figure
        mesh(QQ,TQ,conv_signal.');
        title('Convoluted Signal - Time')
        keyboard
    end
   

end


function [b, a] = swap(a, b)
end

function [num] = roundn(x, n)
    num = n*floor(x/n);
end

function [f] = gaussian(fwhm, x, x0)
    sigma = fwhm/(2*sqrt(2*log(2))); % convert fwhm to std for gaussian pulse
    height = 1/(sigma*sqrt(2*pi)); % gaussian height/ pre factor 
    f = height * exp((-1 .* (x - x0).^2)./ (2*sigma^2));
    f = f./sum(f); % normalise to unity
end 
