classdef gpa < handle
% Geometric Phase Analysis (GPA) Tool in MATLAB
% Written by Subin Bang, MSE SNU @ 2017, quantized@outlook.com
% Any unauthorized use or modification is NOT allowed.

    properties
        
        row;
        col;
        
        input;
        cal_array;
        temp_array;
        
        phases;
        fft_result;
        ifft_result;
        
        PowerSpectrum;
        
        mask_radius;
        mask;
        sigma;
        
        masked_fft;
        raw_phase;
        phase;
        normed_phase;
        ph;

        %angle; % 함수 명과 이름이 겹침
        gx;
        gy;
        dx;
        dy;
        
    end
    
    methods
        function self = gpa(filename)
            
            input_data           = imread(filename);
            self.input           = rgb2gray(input_data);
            self.cal_array       = double(self.input);
            [self.row, self.col] = size(self.input);
            
            self.init_fft();
        end
        
        function result = init_fft(self)
            
            self.temp_array  = self.cal_array;
            self.fft_result  = zeros(self.row, self.col);
            self.ifft_result = zeros(self.row, self.col);
            
            temp            = preFFTshift(self.temp_array);
            self.fft_result = doFFT(temp);
            
            result = self.fft_result;
            self.minrad();
        end
        
        function result = minrad(self)
            % hannwindow = 선택하게 만들기??
            self.temp_array      = self.cal_array;
            hanned               = hannwindow(self.temp_array);
            temp                 = preFFTshift(hanned);
            windowed_fft         = doFFT(temp);
            self.PowerSpectrum   = log10(1+abs(windowed_fft));
            
            x0 = floor(self.col/2);
            y0 = floor(self.row/2);
            
            maxval = real(self.PowerSpectrum(y0,x0));
            
            tic = 1;
            averages = [];
            for r=1:( min(self.row,self.col)/4 )
                vals = [];
                toc = 1;
                
                for i=1:(y0+r)
                    for j=1:(x0+r)
                        dist = distance(x0,y0,j,i);
                        if (dist < r+1 && dist > r-1)
                            vals(toc) = real(self.PowerSpectrum(i,j));
                            toc = toc + 1;
                        end
                    end
                end
                
                [vals,~] = sort(vals,'descend');
                number   = min(numel(vals),10);
                avg      = sum(vals(1:number));
                
                avg           = avg/(number*maxval);
                averages(tic) = avg;
                tic           = tic + 1;
            end
            
            if (numel(averages) < 1)
                result = min(self.row,self.col)/2;
            else % Normalization
                [minval,~] = min(averages);
                averages = averages - minval;
            
                [maxval,~] = max(averages);
                averages = averages/maxval;
            
                inner = 0;
                for i=1:numel(averages)
                    if (averages(i) < 0.9*maxval)
                        inner = i;
                        break;
                    end
                    averages(i) = averages(i) * 0;
                end
                
                for i=1:numel(averages)
                    averages(i) = averages(i) * ( numel(averages) - i - 1 );
                end
                
                [~, I] = max(averages(1+inner:end));
                radius = I(1)+inner;
                result = floor(radius);
                
                %imagesc(abs(self.fft_result));
                self.mask_radius = radius/2;
                self.sigma       = radius/12; % Proper way to estimate sigma?
                self.generate_mask();
            end
            % manual selection for mask radius is needed.
        end
        
        function result = generate_mask(self)
            
            %self.fft_result = self.PowerSpectrum;
            self.temp_array = self.fft_result;
            disp('Select a g-point you want to use for GPA.')
            figure(1)
            %imagesc(abs(self.temp_array));
            imagesc(self.PowerSpectrum);
            [self.gx.first, self.gy.first] = ginput(1);
            close(1)
            
            self.gx.first  = floor(self.gx.first)-self.col/2;
            self.gy.first  = -floor(self.gy.first)+self.row/2;
            x              = (1:self.col)-self.col/2;
            y              = (1:self.row)*(-1)+self.row/2;
            [X, Y]         = meshgrid(x,y);
            
            self.mask.first  = exp( -0.5 * ((X-self.gx.first).^2 + (Y-self.gy.first).^2) / (self.sigma^2) );
            
            disp('Select another g-point you want to use for GPA.')
            figure(1)
            %imagesc(abs(self.temp_array));
            imagesc(self.PowerSpectrum);
            [self.gx.second, self.gy.second] = ginput(1);
            close(1)
            
            self.gx.second  = floor(self.gx.second)-self.col/2;
            self.gy.second  = -floor(self.gy.second)+self.row/2;
            
            self.mask.second = exp( -0.5 * ((X-self.gx.second).^2 + (Y-self.gy.second).^2) / (self.sigma^2) );
            result           = self.mask;
            self.getmaskedfft();
        end
        
        function result = getmaskedfft(self)
            
            self.temp_array        = self.fft_result;
            self.masked_fft.first  = self.temp_array .* self.mask.first;
            self.masked_fft.second = self.temp_array .* self.mask.second;
            
            result = self.masked_fft;
            self.generate_raw_phase();
        end
        
        function result = generate_raw_phase(self)
            
            self.raw_phase.first = doIFFT(self.masked_fft.first);
            self.raw_phase.first = preFFTshift(self.raw_phase.first);
            
            self.raw_phase.second = doIFFT(self.masked_fft.second);
            self.raw_phase.second = preFFTshift(self.raw_phase.second);
            
            self.raw_phase.first  = imag(log(self.raw_phase.first));
            self.raw_phase.second = imag(log(self.raw_phase.second));
            
            result = self.raw_phase;
            self.get_phase();
        end
        
        function result = get_phase(self)
            
            x     = (1:self.col)-self.col/2;
            y     = (1:self.row)*(-1)+self.row/2;
%             x      = 1:self.col;
%             y      = 1:self.row;
            [X,Y] = meshgrid(x,y);
            
%             temp_gx.first  = self.gx.first - self.col/2;
%             temp_gy.first  = -self.gy.first + self.row/2;
%             temp_gx.second = self.gx.second - self.col/2;
%             temp_gy.second = -self.gy.second + self.row/2;
            temp_gx.first = self.gx.first;
            temp_gy.first = self.gy.first;
            temp_gx.second = self.gx.second;
            temp_gy.second = self.gy.second;
            
            self.phase.first  = self.raw_phase.first  - 2*pi* (X*temp_gx.first /self.col + Y*temp_gy.first /self.row);
            self.phase.second = self.raw_phase.second - 2*pi* (X*temp_gx.second/self.col + Y*temp_gy.second/self.row);
            % why divide by self.col and self.row...? = reciprocal lattice vector
            
%             for i=1:self.row
%                 for j=1:self.col
%                     self.phase.first  = self.raw_phase.first  - 2*pi* (j*self.gx.first/self.col  + i*self.gy.first/self.row );
%                     self.phase.second = self.raw_phase.second - 2*pi* (j*self.gx.second/self.col + i*self.gy.second/self.row);
%                     
%                 end
%             end

            result = self.phase;
            self.get_normed_phase();
        end
        
        function result = get_normed_phase(self)
            
            self.normed_phase.first  = self.phase.first  - round( self.phase.first/(2*pi)  ) * (2*pi);
            self.normed_phase.second = self.phase.second - round( self.phase.second/(2*pi) ) * (2*pi);
            
            result = self.normed_phase;
            self.init_distortion();
        end
        
        function result = init_distortion(self)
            
            expPhase.first   = zeros(self.row+2,self.col+2);
            expPhase.second  = zeros(self.row+2,self.col+2);
            dx_kernel.first  = zeros(self.row+2,self.col+2);
            dy_kernel.first  = zeros(self.row+2,self.col+2);
            dx_kernel.second = zeros(self.row+2,self.col+2);
            dy_kernel.second = zeros(self.row+2,self.col+2);
            
            % Fill kernels with pre fft shifted data
            for i=1:3
                dx_kernel.first(i,1)  = 1;
                dx_kernel.first(i,3)  = -1;
                dx_kernel.second(i,1) = 1;
                dx_kernel.second(i,3) = -1;
                
                dy_kernel.first(1,i)  = 1;
                dy_kernel.first(3,i)  = -1;
                dy_kernel.second(1,i) = 1;
                dy_kernel.second(3,i) = -1;
            end
            
            expPhase.first (1:self.row,1:self.col) = exp(self.normed_phase.first  * 1j);
            expPhase.second(1:self.row,1:self.col) = exp(self.normed_phase.second * 1j);
            
            phaseTemp.first  = doFFT(expPhase.first);
            phaseTemp.second = doFFT(expPhase.second);
            xTemp.first      = doFFT(dx_kernel.first);
            yTemp.first      = doFFT(dy_kernel.first);
            xTemp.second     = doFFT(dx_kernel.second);
            yTemp.second     = doFFT(dy_kernel.second);
            
            xTemp.first  = xTemp.first .* phaseTemp.first;
            yTemp.first  = yTemp.first .* phaseTemp.first;
            xTemp.second = xTemp.second .* phaseTemp.second;
            yTemp.second = yTemp.second .* phaseTemp.second;
            
            dx_kernel.first  = doIFFT(xTemp.first);
            dy_kernel.first  = doIFFT(yTemp.first);
            dx_kernel.second = doIFFT(xTemp.second);
            dy_kernel.second = doIFFT(yTemp.second);
            
            self.ph.first  = conj (expPhase.first(2:end-1,2:end-1));
            self.ph.second = conj(expPhase.second(2:end-1,2:end-1));
            
            self.dx.first  = imag( self.ph.first  .* dx_kernel.first (2:end-1,2:end-1) / ( numel(dx_kernel)*6 ) );
            self.dy.first  = imag( self.ph.first  .* dy_kernel.first (2:end-1,2:end-1) / ( numel(dy_kernel)*6 ) );
            self.dx.second = imag( self.ph.second .* dx_kernel.second(2:end-1,2:end-1) / ( numel(dx_kernel)*6 ) );
            self.dy.second = imag( self.ph.second .* dy_kernel.second(2:end-1,2:end-1) / ( numel(dy_kernel)*6 ) );
            
            %Rotation matrix 쓸 것
            
            self.calc_distortion();
        end

%         function result = init_distortion(self)
%             
%         end
        
        function result = calc_distortion(self)
            
            A = [self.gx.first, self.gx.second ; self.gy.first, self.gy.second];
            A = inv(transpose(A));
            %Rotation matrix 곱할 것
            
            factor = -1.0/(2*pi);
            exx = factor * real( A(1,1) * self.dx.first + A(1,2) * self.dx.second );
            exy = factor * real( A(1,1) * self.dy.first + A(1,2) * self.dy.second );
            eyx = factor * real( A(2,1) * self.dx.first + A(2,2) * self.dx.second );
            eyy = factor * real( A(2,1) * self.dy.first + A(2,2) * self.dy.second );
            
            figure(1)
            imagesc(exx);
            figure(2)
            imagesc(exy);
            figure(3)
            imagesc(eyx);
            figure(4)
            imagesc(eyy);
        end
    end
end




function result = doFFT(image, i, j)

if (nargin == 3)
    result = fft2(image,i,j);
elseif (nargin ==2)
    disp('Error occurred during FFT');
else
    result = fft2(image);
end
end

function result = doIFFT(image, i, j)
if (nargin == 3)
    result = ifft2(image,i,j);
elseif (nargin ==2)
    disp('Error occurred during IFFT');
else
    result = ifft2(image);
end
end

function result = preFFTshift(image)

[row,col] = size(image);
result = zeros(row,col);

for i=1:row
    for j=1:col
        result(i,j) = (-1)^(i+j) * image(i,j);
    end
end

end

function result = hannwindow(image)

[row, col] = size(image);
result = zeros(row, col);

arr_y = zeros(row);
arr_x = zeros(col);

for i=1:row
    arr_y(i) = 0.5 * ( 1-cos( (2*pi*i)/(row-1) ) );
end
for i=1:col
    arr_x(i) = 0.5 * ( 1-cos( (2*pi*i)/(col-1) ) );
end

for i=1:row
    for j=1:col
        result(i,j) = image(i,j) * arr_y(i) * arr_x(i);
    end
end
end

function result = distance(x1,y1,x2,y2)
result = sqrt( (x1-x2)^2 + (y1-y2)^2 );
end

    