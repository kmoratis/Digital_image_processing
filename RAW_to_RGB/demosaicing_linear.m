function Ccam = demosaicing_linear(rawim, bayertype, M, N)
% Convert a Bayer-type raw image to RGB, using bilinear interpolation
%  
% Inputs:
%   rawim - the input raw image in Bayer pattern
%   bayertype - a string specifying the Bayer pattern ('RGGB', 'BGGR', 'GBRG', 'GRBG')
%   M, N - the size of the output R, G, B matrices
%
% Output:
%   Ccam - the output matrice containing for R, G, B channels.

    % Determine size of raw image
    [M0, N0, ~] = size(rawim);

    %Determine the ratio of the old dimensions compared to the new dimensions
    y_scale = M/M0;
    x_scale = N/N0;

    % Split raw img to R, G, B channels
    R = NaN(M0,N0);
    G = NaN(M0,N0);
    B = NaN(M0,N0);

    % Create matrices to hold the interpolated R, G, B values
    Rmat = zeros(M,N);
    Gmat = zeros(M,N);
    Bmat = zeros(M,N);

    % Determine the positions of the R, G, B pixels based on the Bayer pattern
    if strcmp(bayertype, 'RGGB') || strcmp(bayertype, 'BGGR')
        % RGGB (BGGR is the same with R, B reversed)
        Rx = 1:2:M0; % only here x is for vertical.. 
        Ry = 1:2:N0;
        Gx1 = 2:2:M0-1; 
        Gy1 = 1:2:N0;
        Gx2 = 1:2:M0;
        Gy2 = 2:2:N0-1;
        Bx = 2:2:M0-1; 
        By = 2:2:N0-1;

        R(Rx,Ry) = rawim(Rx,Ry);
        G(Gx1,Gy1) = rawim(Gx1,Gy1);
        G(Gx2,Gy2) = rawim(Gx2,Gy2);
        B(Bx,By) = rawim(Bx,By);

        %% R mat ( R channel )
        % Fill in the R mat with the existing values 
        for count1 = 0:2:M-1
            for count2 = 0:2:N-1
                y = min(max(1+round(count1./y_scale),1),M0); % vertical
                x = min(max(1+round(count2./x_scale),1),N0); % horizontal
 
                if(~isnan(R(y,x)))
                    Rmat(count1+1,count2+1) = R(y,x);
                else
                    Rmat(count1+1,count2+1) = 0;
                end
            end
        end

        % Fill the rest of the R mat (previous G, B pixel spots) using bilinear interpolation
        % previous G1 positions
        % linear interpolation on x axis
        for count1 = 0:2:M-1
            for count2 = 1:2:N-2
                % actual coordinates
                y = min(max(1+count1/y_scale,1),M0); % vertical
                x = min(max(1+count2/x_scale,1),N0); % horizontal
         
                % left neighbor
                y1 = min(max(1+round(count1/y_scale),1),M0); % vertical
                x1 = min(max(1+round(count2/x_scale) - 1,1),N0); % horizontal
                Q1 = R(y1,x1);

                % right neighbor
                y2 = min(max(1+round(count1/y_scale),1),M0); % vertical
                x2 = min(max(1+round(count2/x_scale) + 1,1),N0); % horizontal
                Q2 = R(y2,x2);

                Qint = ((x2-x)/(x2-x1))*Q1 + ((x-x1)/(x2-x1))*Q2;
                Rmat(count1+1, count2+1) = Qint;
            end
        end

        % previous G2 positions
        % linear interpolation on y axis
        for count1 = 1:2:M-2
            for count2 = 0:2:N-1
                % actual coordinates
                y = min(max(1+count1/y_scale,1),M0); % vertical
                x = min(max(1+count2/x_scale,1),N0); % horizontal
         
                % upper neighbor
                y1 = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                x1 = min(max(1+round(count2./x_scale),1),N0); % horizontal
                Q1 = R(y1,x1);

                % lower neighbor
                y2 = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                x2 = min(max(1+round(count2./x_scale),1),N0); % horizontal
                Q2 = R(y2,x2);

                Qint = ((y2-y)/(y2-y1))*Q1 + ((y-y1)/(y2-y1))*Q2;
                Rmat(count1+1, count2+2) = Qint;
            end
        end

        % previous B positions
        % bilinear interpolation
        for count1 = 1:2:M-2
            for count2 = 1:2:N-2
                % actual coordinates
                y = min(max(1+count1/y_scale,1),M0); % vertical
                x = min(max(1+count2/x_scale,1),N0); % horizontal
         
                % coordinates for neighbors
                x1 = min(max(1+round(count2./x_scale) - 1,1),N0);
                x2 = min(max(1+round(count2./x_scale) + 1,1),N0);
                y1 = min(max(1+round(count1./y_scale) - 1,1),M0);
                y2 = min(max(1+round(count1./y_scale) + 1,1),M0);

                Q11 = R(y1,x1); % upper left neighbor
                Q12 = R(y1,x2); % upper right neighbor
                Q21 = R(y2,x1); % lower left neighbor
                Q22 = R(y2,x2); % lower right neighbor

                % interpolating on 2 first axis and the result between them
                Qu = ((x2-x)/(x2-x1))*Q11 + ((x-x1)/(x2-x1))*Q12;
                Ql = ((x2-x)/(x2-x1))*Q21 + ((x-x1)/(x2-x1))*Q22;
                Qint = ((y2-y)/(y2-y1))*Qu + ((y-y1)/(y2-y1))*Ql;
                Rmat(count1+1, count2+2) = Qint;
            end
        end

        %% G mat ( G channel )
        % Fill in the G mat with the existing values
        % G1
        for count1 = 0:2:M-1
            for count2 = 1:2:N-2
                y = min(max(1+round(count1./y_scale),1),M0); % vertical
                x = min(max(1+round(count2./x_scale),1),N0); % horizontal

                if(~isnan(G(y,x)))
                    Gmat(count1+1,count2+1) = G(y,x);
                else
                    Gmat(count1+1,count2+1) = 0;
                end
            end
        end

        % G2
        for count1 = 1:2:M-2
            for count2 = 0:2:N-1
                y = min(max(1+round(count1./y_scale),1),M0); % vertical
                x = min(max(1+round(count2./x_scale),1),N0); % horizontal

                if(~isnan(G(y,x)))
                    Gmat(count1+1,count2+1) = G(y,x);
                else
                    Gmat(count1+1,count2+1) = 0;
                end
            end
        end

        % Fill the rest of the G mat (previous R, B pixel spots) using bilinear interpolation
        % previous R positions
        % linear interpolation on x and y axis and taking the mean of the two
        for count1 = 0:2:M-1
            for count2 = 0:2:N-1
                % actual coordinates
                y = min(max(1+count1/y_scale,1),M0); % vertical
                x = min(max(1+count2/x_scale,1),N0); % horizontal
         
                % coordinates for neighbors
                x1 = min(max(1+round(count2./x_scale) - 1,1),N0);
                x2 = min(max(1+round(count2./x_scale) + 1,1),N0);
                y1 = min(max(1+round(count1./y_scale) - 1,1),M0);
                y2 = min(max(1+round(count1./y_scale) + 1,1),M0);
                y0 = min(max(1+round(count1./y_scale),1),M0);
                x0 = min(max(1+round(count2./x_scale),1),N0);

                Q11 = G(y1,x0);  % upper neighbor
                Q12 = G(y2,x0); % lower neighbor
                Q21 = G(y0,x1); % left neighbor
                Q22 = G(y0,x2); % right neighbor
                
                % interpolating on x axis and y axis and taking the average
                Qx = ((x2-x)/(x2-x1))*Q21 + ((x-x1)/(x2-x1))*Q22;
                Qy = ((y2-y)/(y2-y1))*Q11 + ((y-y1)/(y2-y1))*Q12;

                Gmat(count1+1, count2+1) = mean(nonzeros([Qx Qy]), 'omitnan');
            end
        end

        % previous B positions
        % linear interpolation on x and y axis and taking the mean of the two
        for count1 = 1:2:M-2
            for count2 = 1:2:N-2
                % actual coordinates
                y = min(max(1+count1/y_scale,1),M0); % vertical
                x = min(max(1+count2/x_scale,1),N0); % horizontal
         
                % coordinates for neighbors
                x1 = min(max(1+round(count2./x_scale) - 1,1),N0);
                x2 = min(max(1+round(count2./x_scale) + 1,1),N0);
                y1 = min(max(1+round(count1./y_scale) - 1,1),M0);
                y2 = min(max(1+round(count1./y_scale) + 1,1),M0);
                y0 = min(max(1+round(count1./y_scale),1),M0);
                x0 = min(max(1+round(count2./x_scale),1),N0);

                Q11 = G(y1,x0);  % upper neighbor
                Q12 = G(y2,x0); % lower neighbor
                Q21 = G(y0,x1); % left neighbor
                Q22 = G(y0,x2); % right neighbor
                
                % interpolating on x axis and y axis and taking the average
                Qx = ((x2-x)/(x2-x1))*Q21 + ((x-x1)/(x2-x1))*Q22;
                Qy = ((y2-y)/(y2-y1))*Q11 + ((y-y1)/(y2-y1))*Q12;

                Gmat(count1+1, count2+1) = mean(nonzeros([Qx Qy]), 'omitnan');
            end
        end

        %% B mat ( B channel )
        % Fill in the B mat with the existing values
        for count1 = 1:2:M-2
            for count2 = 1:2:N-2
                y = min(max(1+round(count1./y_scale),1),M0); % vertical
                x = min(max(1+round(count2./x_scale),1),N0); % horizontal
                % if might not be necessary (2)
                if(~isnan(B(y,x)))
                    Bmat(count1+1,count2+1) = B(y,x);
                else
                   Bmat(count1+1,count2+1) = 0;
                end
            end
        end

        % Fill the rest of the B mat (previous R, G pixel spots) using bilinear interpolation
        % previous R positions
        % bilinear interpolation
        for count1 = 0:2:M-1
            for count2 = 0:2:N-1
                % actual coordinates
                y = min(max(1+count1/y_scale,1),M0); % vertical
                x = min(max(1+count2/x_scale,1),N0); % horizontal
         
                % coordinates for neighbors
                x1 = min(max(1+round(count2./x_scale) - 1,1),N0);
                x2 = min(max(1+round(count2./x_scale) + 1,1),N0);
                y1 = min(max(1+round(count1./y_scale) - 1,1),M0);
                y2 = min(max(1+round(count1./y_scale) + 1,1),M0);

                Q11 = B(y1,x1); % upper left neighbor
                Q12 = B(y1,x2); % upper right neighbor
                Q21 = B(y2,x1); % lower left neighbor
                Q22 = B(y2,x2); % lower right neighbor

                % interpolating on 2 first axis and the result between them
                Qu = ((x2-x)/(x2-x1))*Q11 + ((x-x1)/(x2-x1))*Q12;
                Ql = ((x2-x)/(x2-x1))*Q21 + ((x-x1)/(x2-x1))*Q22;
                Qint = ((y2-y)/(y2-y1))*Qu + ((y-y1)/(y2-y1))*Ql;
                Bmat(count1+1, count2+2) = Qint;
            end
        end

        % previous G1 positions
        % linear interpolation on y axis
        for count1 = 0:2:M-1
            for count2 = 1:2:N-2
                % actual coordinates
                y = min(max(1+count1/y_scale,1),M0); % vertical
                x = min(max(1+count2/x_scale,1),N0); % horizontal
         
                % upper neighbor
                y1 = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                x1 = min(max(1+round(count2./x_scale),1),N0); % horizontal
                Q1 = B(y1,x1);

                % lower neighbor
                y2 = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                x2 = min(max(1+round(count2./x_scale),1),N0); % horizontal
                Q2 = B(y2,x2);

                Qint = ((y2-y)/(y2-y1))*Q1 + ((y-y1)/(y2-y1))*Q2;
                Bmat(count1+1, count2+2) = Qint;
            end
        end

        % previous G2 positions
        % linear interpolation on x axis
        for count1 = 1:2:M-2
            for count2 = 0:2:N-1
                % actual coordinates
                y = min(max(1+count1/y_scale,1),M0); % vertical
                x = min(max(1+count2/x_scale,1),N0); % horizontal
         
                % left neighbor
                y1 = min(max(1+round(count1./y_scale),1),M0); % vertical
                x1 = min(max(1+round(count2./x_scale) - 1,1),N0); % horizontal
                Q1 = B(y1,x1);

                % right neighbor
                y2 = min(max(1+round(count1./y_scale),1),M0); % vertical
                x2 = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                Q2 = B(y2,x2);

                Qint = ((x2-x)/(x2-x1))*Q1 + ((x-x1)/(x2-x1))*Q2;
                Bmat(count1+1, count2+1) = Qint;
            end
        end
    
        % dealing with resizing (works fine for M0*N0)
        if (N ~= N0) % Fill columns that don't have values
            for count2 = 0:1:N-2
                % Red
                if (~any(Rmat(1,count2+1)) || ~any(Rmat(2,count2+1)) || anynan(Rmat(1, count2+1)) || anynan(Rmat(2, count2+1))) % if first 2 elements of a column are zeros then replace it with interpolation
                    for count1 = 0:1:M-1
                        if count2 == 0
                            neighbors = [Rmat(count1+1, 1), Rmat(count1+1, 2)];
                        else
                            neighbors = [Rmat(count1+1, count2), Rmat(count1+1, count2+2)];
                        end
                        Rmat(count1+1, count2+1) = mean(nonzeros(neighbors), 'omitnan');
                    end
                end
                % Green
                if (~any(Gmat(1,count2+1)) || ~any(Gmat(2,count2+1)) || anynan(Gmat(1, count2+1)) || anynan(Gmat(2, count2+1))) % if first 2 elements of a column are zeros then replace it with interpolation
                    for count1 = 0:1:M-1
                        if count2 == 0
                            neighbors = [Gmat(count1+1, 1), Gmat(count1+1, 2)];
                        else
                            neighbors = [Gmat(count1+1, count2), Gmat(count1+1, count2+2)];
                        end
                        Gmat(count1+1, count2+1) = mean(nonzeros(neighbors), 'omitnan');
                    end
                end
                % Blue
                if (~any(Bmat(1,count2+1)) || ~any(Bmat(2,count2+1)) || anynan(Bmat(1, count2+1)) || anynan(Bmat(2, count2+1))) % if first 2 elements of a column are zeros then replace it with interpolation
                    for count1 = 0:1:M-1
                        if count2 == 0
                            neighbors = [Bmat(count1+1, 1), Bmat(count1+1, 2)];
                        else
                            neighbors = [Bmat(count1+1, count2), Bmat(count1+1, count2+2)];
                        end
                        Bmat(count1+1, count2+1) = mean(nonzeros(neighbors), 'omitnan');
                    end
                end
            end
        end

        if (M ~= M0) % Fill rows that don't have values
            for count1 = 0:1:M-2
                % Red
                if (~any(Rmat(count1+1, 1)) || ~any(Rmat(count1+1, 2)) || anynan(Rmat(count1+1, 1)) || anynan(Rmat(count1+1, 2))) % if first 2 elements of a row are zeros then replace it with interpolation
                    for count2 = 0:1:N-1
                        if count1 == 0
                            neighbors = [Rmat(1, count2+1), Rmat(2, count2+1)];
                        else
                            neighbors = [Rmat(count1, count2+1), Rmat(count1+2, count2+1)];
                        end
                        Rmat(count1+1, count2+1) = mean(nonzeros(neighbors), 'omitnan');
                    end
                end
                % Green
                if (~any(Gmat(count1+1, 2)) || ~any(Gmat(count1+1, 3)) || anynan(Gmat(count1+1, 2))) || anynan(Gmat(count1+1, 2)) % if elements 2,3  of a row are zeros then replace it with interpolation
                    for count2 = 0:1:N-1
                        if count1 == 0
                            neighbors = [Gmat(1, count2+1), Gmat(2, count2+1)];
                        else
                            neighbors = [Gmat(count1, count2+1), Gmat(count1+2, count2+1)];
                        end
                        Gmat(count1+1, count2+1) = mean(nonzeros(neighbors), 'omitnan');
                    end
                end
                % Blue
                if (~any(Bmat(count1+1, 2)) || ~any(Bmat(count1+1, 3)) || anynan(Bmat(count1+1, 2))) || anynan(Bmat(count1+1, 2)) % if elements 2,3  of a row are zeros then replace it with interpolation
                    for count2 = 0:1:N-1
                        if count1 == 0
                            neighbors = [Bmat(1, count2+1), Bmat(2, count2+1)];
                        else
                            neighbors = [Bmat(count1, count2+1), Bmat(count1+2, count2+1)];
                        end
                        Bmat(count1+1, count2+1) = mean(nonzeros(neighbors), 'omitnan');
                    end
                end
            end
        end
        

        % image limits handling
        Rmat(M,:) = Rmat(M-1,:);
        Rmat(2:2:M,1) = Rmat(1:2:M,1);
        Rmat(1:2:M,N) = Rmat(2:2:M,N);

        Gmat(M,:) = Gmat(M-1,:);
        Gmat(:,N) = Gmat(:,N-1);

        Bmat(1,:) = Bmat(2,:);
        Bmat(M,:) = Bmat(M-2,:);
        Bmat(M-1,:) = Bmat(M-2,:);
        Bmat(1:2:M,2) = Bmat(1:2:M,3);
        Bmat(:,1) = Bmat(:,2);
        Bmat(2:2:M,N-1) = Bmat(2:2:M,N-2);
        Bmat(:,N) = Bmat(:,N-1);

        if strcmp(bayertype, 'RGGB')
            %Rmat = zeros(M,N);
            %Gmat = zeros(M,N);
            Ccam = cat(3, Rmat, Gmat, Bmat);
        else % if 'BGGR'
            Ccam = cat(3, Bmat, Gmat, Rmat);
        end
    
    elseif strcmp(bayertype, 'GBRG') || strcmp(bayertype, 'GRBG')
        % Determine the positions of the R, G, B pixels based on the Bayer pattern
        % GRBG (GBRG is the same with R,B reversed)
        Rx = 1:2:M0; % only here x is for vertical.. 
        Ry = 2:2:N0-1;
        Gx1 = 1:2:M0; 
        Gy1 = 1:2:N0;
        Gx2 = 2:2:M0-1;
        Gy2 = 2:2:N0-1;
        Bx = 2:2:M0-1; 
        By = 1:2:N0;

        R(Rx,Ry) = rawim(Rx,Ry);
        G(Gx1,Gy1) = rawim(Gx1,Gy1);
        G(Gx2,Gy2) = rawim(Gx2,Gy2);
        B(Bx,By) = rawim(Bx,By);

        %% R mat ( R channel )
        % Fill in the R mat with the existing values 
        for count1 = 0:2:M-1
            for count2 = 1:2:N-2
                y = min(max(1+round(count1./y_scale),1),M0); % vertical
                x = min(max(1+round(count2./x_scale),1),N0); % horizontal
 
                if(~isnan(R(y,x)))
                    Rmat(count1+1,count2+1) = R(y,x);
                else
                    Rmat(count1+1,count2+1) = 0;
                end
            end
        end

        % Fill the rest of the R mat (previous G, B pixel spots) using bilinear interpolation
        % previous G1 positions
        % linear interpolation on x axis
        for count1 = 0:2:M-1
            for count2 = 0:2:N-1
                % actual coordinates
                y = min(max(1+count1/y_scale,1),M0); % vertical
                x = min(max(1+count2/x_scale,1),N0); % horizontal
         
                % left neighbor
                y1 = min(max(1+round(count1./y_scale),1),M0); % vertical
                x1 = min(max(1+round(count2./x_scale) - 1,1),N0); % horizontal
                Q1 = R(y1,x1);

                % right neighbor
                y2 = min(max(1+round(count1./y_scale),1),M0); % vertical
                x2 = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                Q2 = R(y2,x2);

                Qint = ((x2-x)/(x2-x1))*Q1 + ((x-x1)/(x2-x1))*Q2;
                Rmat(count1+1, count2+1) = Qint;
            end
        end

        % previous G2 positions
        % linear interpolation on y axis
        for count1 = 1:2:M-2
            for count2 = 1:2:N-2
                % actual coordinates
                y = min(max(1+count1/y_scale,1),M0); % vertical
                x = min(max(1+count2/x_scale,1),N0); % horizontal
         
                % upper neighbor
                y1 = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                x1 = min(max(1+round(count2./x_scale),1),N0); % horizontal
                Q1 = R(y1,x1);

                % lower neighbor
                y2 = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                x2 = min(max(1+round(count2./x_scale),1),N0); % horizontal
                Q2 = R(y2,x2);

                Qint = ((y2-y)/(y2-y1))*Q1 + ((y-y1)/(y2-y1))*Q2;
                Rmat(count1+1, count2+2) = Qint;
            end
        end

        % previous B positions
        % bilinear interpolation
        for count1 = 1:2:M-2
            for count2 = 0:2:N-1
                % actual coordinates
                y = min(max(1+count1/y_scale,1),M0); % vertical
                x = min(max(1+count2/x_scale,1),N0); % horizontal
         
                % coordinates for neighbors
                x1 = min(max(1+round(count2./x_scale) - 1,1),N0);
                x2 = min(max(1+round(count2./x_scale) + 1,1),N0);
                y1 = min(max(1+round(count1./y_scale) - 1,1),M0);
                y2 = min(max(1+round(count1./y_scale) + 1,1),M0);

                Q11 = R(y1,x1); % upper left neighbor
                Q12 = R(y1,x2); % upper right neighbor
                Q21 = R(y2,x1); % lower left neighbor
                Q22 = R(y2,x2); % lower right neighbor

                % interpolating on 2 first axis and the result between them
                Qu = ((x2-x)/(x2-x1))*Q11 + ((x-x1)/(x2-x1))*Q12;
                Ql = ((x2-x)/(x2-x1))*Q21 + ((x-x1)/(x2-x1))*Q22;
                Qint = ((y2-y)/(y2-y1))*Qu + ((y-y1)/(y2-y1))*Ql;
                Rmat(count1+1, count2+2) = Qint;
            end
        end

        %% G mat ( G channel )
        % Fill in the G mat with the existing values
        % G1
        for count1 = 0:2:M-1
            for count2 = 0:2:N-1
                y = min(max(1+round(count1./y_scale),1),M0); % vertical
                x = min(max(1+round(count2./x_scale),1),N0); % horizontal

                if(~isnan(G(y,x)))
                    Gmat(count1+1,count2+1) = G(y,x);
                else
                    Gmat(count1+1,count2+1) = 0;
                end
            end
        end

        % G2
        for count1 = 1:2:M-2
            for count2 = 1:2:N-2
                y = min(max(1+round(count1./y_scale),1),M0); % vertical
                x = min(max(1+round(count2./x_scale),1),N0); % horizontal

                if(~isnan(G(y,x)))
                    Gmat(count1+1,count2+1) = G(y,x);
                else
                    Gmat(count1+1,count2+1) = 0;
                end
            end
        end

        % Fill the rest of the G mat (previous R, B pixel spots) using nearest neighbor interpolation
        % previous R positions
        % interpolation on x and y axis and taking the average of the two
        for count1 = 0:2:M-1
            for count2 = 1:2:N-2
                % actual coordinates
                y = min(max(1+count1/y_scale,1),M0); % vertical
                x = min(max(1+count2/x_scale,1),N0); % horizontal
         
                % coordinates for neighbors
                x1 = min(max(1+round(count2./x_scale) - 1,1),N0);
                x2 = min(max(1+round(count2./x_scale) + 1,1),N0);
                y1 = min(max(1+round(count1./y_scale) - 1,1),M0);
                y2 = min(max(1+round(count1./y_scale) + 1,1),M0);
                y0 = min(max(1+round(count1./y_scale),1),M0);
                x0 = min(max(1+round(count2./x_scale),1),N0);

                Q11 = G(y1,x0);  % upper neighbor
                Q12 = G(y2,x0); % lower neighbor
                Q21 = G(y0,x1); % left neighbor
                Q22 = G(y0,x2); % right neighbor
                
                % interpolating on x axis and y axis and taking the average
                Qx = ((x2-x)/(x2-x1))*Q21 + ((x-x1)/(x2-x1))*Q22;
                Qy = ((y2-y)/(y2-y1))*Q11 + ((y-y1)/(y2-y1))*Q12;

                Gmat(count1+1, count2+1) = mean(nonzeros([Qx Qy]), 'omitnan');
            end
        end

        % previous B positions
        % interpolation on x and y axis and taking the average of the two
        for count1 = 1:2:M-2
            for count2 = 0:2:N-1
                % actual coordinates
                y = min(max(1+count1/y_scale,1),M0); % vertical
                x = min(max(1+count2/x_scale,1),N0); % horizontal
         
                % coordinates for neighbors
                x1 = min(max(1+round(count2./x_scale) - 1,1),N0);
                x2 = min(max(1+round(count2./x_scale) + 1,1),N0);
                y1 = min(max(1+round(count1./y_scale) - 1,1),M0);
                y2 = min(max(1+round(count1./y_scale) + 1,1),M0);
                y0 = min(max(1+round(count1./y_scale),1),M0);
                x0 = min(max(1+round(count2./x_scale),1),N0);

                Q11 = G(y1,x0);  % upper neighbor
                Q12 = G(y2,x0); % lower neighbor
                Q21 = G(y0,x1); % left neighbor
                Q22 = G(y0,x2); % right neighbor
                
                % interpolating on x axis and y axis and taking the average
                Qx = ((x2-x)/(x2-x1))*Q21 + ((x-x1)/(x2-x1))*Q22;
                Qy = ((y2-y)/(y2-y1))*Q11 + ((y-y1)/(y2-y1))*Q12;

                Gmat(count1+1, count2+1) = mean(nonzeros([Qx Qy]), 'omitnan');
            end
        end

        %% B mat ( B channel )
        % Fill in the B mat with the existing values
        for count1 = 1:2:M-2
            for count2 = 0:2:N-1
                y = min(max(1+round(count1./y_scale),1),M0); % vertical
                x = min(max(1+round(count2./x_scale),1),N0); % horizontal
                % if might not be necessary (2)
                if(~isnan(B(y,x)))
                    Bmat(count1+1,count2+1) = B(y,x);
                else
                   Bmat(count1+1,count2+1) = 0;
                end
            end
        end

        % Fill the rest of the B mat (previous R, G pixel spots) using nearest neighbor interpolation
        % previous R positions
        % bilinear interpolation
        for count1 = 0:2:M-1
            for count2 = 1:2:N-2
                % actual coordinates
                y = min(max(1+count1/y_scale,1),M0); % vertical
                x = min(max(1+count2/x_scale,1),N0); % horizontal
         
                % coordinates for neighbors
                x1 = min(max(1+round(count2./x_scale) - 1,1),N0);
                x2 = min(max(1+round(count2./x_scale) + 1,1),N0);
                y1 = min(max(1+round(count1./y_scale) - 1,1),M0);
                y2 = min(max(1+round(count1./y_scale) + 1,1),M0);

                Q11 = B(y1,x1); % upper left neighbor
                Q12 = B(y1,x2); % upper right neighbor
                Q21 = B(y2,x1); % lower left neighbor
                Q22 = B(y2,x2); % lower right neighbor

                % interpolating on 2 first axis and the result between them
                Qu = ((x2-x)/(x2-x1))*Q11 + ((x-x1)/(x2-x1))*Q12;
                Ql = ((x2-x)/(x2-x1))*Q21 + ((x-x1)/(x2-x1))*Q22;
                Qint = ((y2-y)/(y2-y1))*Qu + ((y-y1)/(y2-y1))*Ql;
                Bmat(count1+1, count2+2) = Qint;
            end
        end

        % previous G1 positions
        % linear interpolation on y axis
        for count1 = 0:2:M-1
            for count2 = 0:2:N-1
                % actual coordinates
                y = min(max(1+count1/y_scale,1),M0); % vertical
                x = min(max(1+count2/x_scale,1),N0); % horizontal
         
                % upper neighbor
                y1 = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                x1 = min(max(1+round(count2./x_scale),1),N0); % horizontal
                Q1 = B(y1,x1);

                % lower neighbor
                y2 = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                x2 = min(max(1+round(count2./x_scale),1),N0); % horizontal
                Q2 = B(y2,x2);

                Qint = ((y2-y)/(y2-y1))*Q1 + ((y-y1)/(y2-y1))*Q2;
                Bmat(count1+1, count2+2) = Qint;
            end
        end

        % previous G2 positions
        % linear interpolation on x axis
        for count1 = 1:2:M-2
            for count2 = 1:2:N-2
                % actual coordinates
                y = min(max(1+count1/y_scale,1),M0); % vertical
                x = min(max(1+count2/x_scale,1),N0); % horizontal
         
                % left neighbor
                y1 = min(max(1+round(count1./y_scale),1),M0); % vertical
                x1 = min(max(1+round(count2./x_scale) - 1,1),N0); % horizontal
                Q1 = B(y1,x1);

                % right neighbor
                y2 = min(max(1+round(count1./y_scale),1),M0); % vertical
                x2 = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                Q2 = B(y2,x2);

                Qint = ((x2-x)/(x2-x1))*Q1 + ((x-x1)/(x2-x1))*Q2;
                Bmat(count1+1, count2+1) = Qint;
            end
        end

        % dealing with resizing (works fine for M0*N0)
        if (N ~= N0) % Fill columns that don't have values
            for count2 = 0:1:N-2
                % Red
                if (~any(Rmat(1,count2+1)) || ~any(Rmat(2,count2+1)) || anynan(Rmat(2:3, count2+1))) % if first 2 elements of a column are zeros then replace it with interpolation
                    for count1 = 0:1:M-1
                        if count2 == 0
                            neighbors = [Rmat(count1+1, 1), Rmat(count1+1, 2)];
                        else
                            neighbors = [Rmat(count1+1, count2), Rmat(count1+1, count2+2)];
                        end
                        Rmat(count1+1, count2+1) = mean(nonzeros(neighbors), 'omitnan');
                    end
                end
                % Green
                if (~any(Gmat(1,count2+1)) || ~any(Gmat(2,count2+1)) || anynan(Gmat(2:3, count2+1))) % if first 2 elements of a column are zeros then replace it with interpolation
                    for count1 = 0:1:M-1
                        if count2 == 0
                            neighbors = [Gmat(count1+1, 1), Gmat(count1+1, 2)];
                        else
                            neighbors = [Gmat(count1+1, count2), Gmat(count1+1, count2+2)];
                        end
                        Gmat(count1+1, count2+1) = mean(nonzeros(neighbors), 'omitnan');
                    end
                end
                % Blue
                if (~any(Bmat(1,count2+1)) || ~any(Bmat(2,count2+1)) || anynan(Bmat(2:3, count2+1))) % if first 2 elements of a column are zeros then replace it with interpolation
                    for count1 = 0:1:M-1
                        if count2 == 0
                            neighbors = [Bmat(count1+1, 1), Bmat(count1+1, 2)];
                        else
                            neighbors = [Bmat(count1+1, count2), Bmat(count1+1, count2+2)];
                        end
                        Bmat(count1+1, count2+1) = mean(nonzeros(neighbors), 'omitnan');
                    end
                end
            end
        end

        if (M ~= M0) % Fill rows that don't have values
            for count1 = 0:1:M-2
                % Red
                if (~any(Rmat(count1+1, 1)) || ~any(Rmat(count1+1, 2)) || anynan(Rmat(count1+1, 1)) || anynan(Rmat(count1+1, 2))) % if first 2 elements of a row are zeros then replace it with interpolation
                    for count2 = 0:1:N-1
                        if count1 == 0
                            neighbors = [Rmat(1, count2+1), Rmat(2, count2+1)];
                        else
                            neighbors = [Rmat(count1, count2+1), Rmat(count1+2, count2+1)];
                        end
                        Rmat(count1+1, count2+1) = mean(nonzeros(neighbors), 'omitnan');
                    end
                end
                % Green
                if (~any(Gmat(count1+1, 2)) || ~any(Gmat(count1+1, 3)) || anynan(Gmat(count1+1, 2))) || anynan(Gmat(count1+1, 2)) % if elements 2,3  of a row are zeros then replace it with interpolation
                    for count2 = 0:1:N-1
                        if count1 == 0
                            neighbors = [Gmat(1, count2+1), Gmat(2, count2+1)];
                        else
                            neighbors = [Gmat(count1, count2+1), Gmat(count1+2, count2+1)];
                        end
                        Gmat(count1+1, count2+1) = mean(nonzeros(neighbors), 'omitnan');
                    end
                end
                % Blue
                if (~any(Bmat(count1+1, 2)) || ~any(Bmat(count1+1, 3)) || anynan(Bmat(count1+1, 2))) || anynan(Bmat(count1+1, 2)) % if elements 2,3  of a row are zeros then replace it with interpolation
                    for count2 = 0:1:N-1
                        if count1 == 0
                            neighbors = [Bmat(1, count2+1), Bmat(2, count2+1)];
                        else
                            neighbors = [Bmat(count1, count2+1), Bmat(count1+2, count2+1)];
                        end
                        Bmat(count1+1, count2+1) = mean(nonzeros(neighbors), 'omitnan');
                    end
                end
            end
        end

        % image limits handling
        Rmat(M,:) = Rmat(M-1,:);
        Rmat(:,N) = Rmat(:,N-1);

        Gmat(:,1) = Gmat(:,2);
        Gmat(M,:) = Gmat(M-1,:);
        Gmat(:,N) = Gmat(:,N-1);

        Bmat(M,:) = Bmat(M-1,:);
        Bmat(:,N) = Bmat(:,N-1);
    
        if strcmp(bayertype, 'GRBG')
            %Bmat = zeros(M,N);
            %Rmat = zeros(M,N);
            Ccam = cat(3, Rmat, Gmat, Bmat);
        else % if 'GBRG'
            Ccam = cat(3, Bmat, Gmat, Rmat);
        end

    else
        error('Invalid Bayer pattern type');
    end
end
