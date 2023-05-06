function Ccam = demosaicing_nearest(rawim, bayertype, M, N)
% Convert a Bayer-type raw image to RGB, using nearest neighbor interpolation
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
    
        % Split raw img to R, G, B channels
        R = NaN(M0,N0);
        G = NaN(M0,N0);
        B = NaN(M0,N0);

        R(Rx,Ry) = rawim(Rx,Ry);
        G(Gx1,Gy1) = rawim(Gx1,Gy1);
        G(Gx2,Gy2) = rawim(Gx2,Gy2);
        B(Bx,By) = rawim(Bx,By);

        % Create matrices to hold the interpolated R, G, B values
        Rmat = zeros(M,N);
        Gmat = zeros(M,N);
        Bmat = zeros(M,N);

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

        % Fill the rest of the R mat (previous G, B pixel spots) using nearest neighbor interpolation
        % previous G1 positions
        for count1 = 0:2:M-1
            for count2 = 1:2:N-2
                % before round value was closer to the left neighbour
                if(round(count2./x_scale) > count2/x_scale)
                    y = min(max(1+round(count1./y_scale),1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) - 1,1),N0); % horizontal
                else
                    y = min(max(1+round(count1./y_scale),1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                end

                if(~isnan(R(y,x)))
                    Rmat(count1+1,count2+1) = R(y,x);
                else
                    Rmat(count1+1,count2+1) = 0;
                end
            end
        end

        % previous G2 positions
        for count1 = 1:2:M-2
            for count2 = 0:2:N-1
                % before round value was closer to the upper neighbour
                if(round(count1/y_scale) > count1/y_scale)
                    y = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale),1),N0); % horizontal
                else
                    y = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale),1),N0); % horizontal
                end

                if(~isnan(R(y,x)))
                    Rmat(count1+1,count2+1) = R(y,x);
                else
                    Rmat(count1+1,count2+1) = 0;
                end
            end
        end

        % previous B positions
        for count1 = 1:2:M-2
            for count2 = 1:2:N-2
                % before round value was closer to the left and upper neighbour
                if(round(count2/x_scale) > count2/x_scale && round(count1/y_scale) > count1/y_scale)
                    y = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) -1,1),N0); % horizontal
                % before round value was closer to the right and upper neighbour
                elseif(round(count2/x_scale) < count2/x_scale && round(count1/y_scale) > count1/y_scale)
                    y = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                % before round value was closer to the left and lower neighbour
                elseif(round(count2/x_scale) > count2/x_scale && round(count1/y_scale) < count1/y_scale)
                    y = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) - 1,1),N0); % horizontal
                % before round value was closer to the right and lower neighbour
                else
                    y = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                end

                if(~isnan(R(y,x)))
                    Rmat(count1+1,count2+1) = R(y,x);
                else
                    Rmat(count1+1,count2+1) = 0;
                end
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

        % Fill the rest of the G mat (previous R, B pixel spots) using nearest neighbor interpolation
        % previous R positions
        for count1 = 0:2:M-1
            for count2 = 0:2:N-1
                % distances from borders: left, right, upper, lower 
                dist = [abs(round(count2/x_scale)-count2/x_scale), 1 - abs(round(count2/x_scale)-count2/x_scale), abs(round(count1/y_scale)-count1/y_scale), 1 - abs(round(count1/y_scale)-count1/y_scale)];
                [~,I] = min(dist);
                if(I == 1) % left is the nearest
                    y = min(max(1+round(count1./y_scale),1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) - 1,1),N0); % horizontal
                elseif(I == 2) % right
                    y = min(max(1+round(count1./y_scale),1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                elseif(I == 3) % upper
                    y = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale),1),N0); % horizontal
                else % lower
                    y = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale),1),N0); % horizontal
                end

                if(~isnan(G(y,x)))
                    Gmat(count1+1,count2+1) = G(y,x);
                else
                    Gmat(count1+1,count2+1) = 0;
                end
            end
        end

        % previous B positions
        for count1 = 1:2:M-2
            for count2 = 1:2:N-2
                % distances from borders: left, right, upper, lower 
                dist = [abs(round(count2/x_scale)-count2/x_scale), 1 - abs(round(count2/x_scale)-count2/x_scale), abs(round(count1/y_scale)-count1/y_scale), 1 - abs(round(count1/y_scale)-count1/y_scale)];
                [~,I] = min(dist);
                if(I == 1) % left is the nearest
                    y = min(max(1+round(count1./y_scale),1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) - 1,1),N0); % horizontal
                elseif(I == 2) % right
                    y = min(max(1+round(count1./y_scale),1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                elseif(I == 3) % upper
                    y = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale),1),N0); % horizontal
                else % lower
                    y = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale),1),N0); % horizontal
                end

                % if might not be necessary (2)
                if(~isnan(G(y,x)))
                    Gmat(count1+1,count2+1) = G(y,x);
                else
                    Gmat(count1+1,count2+1) = 0;
                end
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

        % Fill the rest of the B mat (previous R, G pixel spots) using nearest neighbor interpolation
        % previous R positions
        for count1 = 0:2:M-1
            for count2 = 0:2:N-1
                % before round value was closer to the left and upper neighbour
                if(round(count2/x_scale) > count2/x_scale && round(count1/y_scale) > count1/y_scale)
                    y = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) -1,1),N0); % horizontal
                % before round value was closer to the right and upper neighbour
                elseif(round(count2/x_scale) < count2/x_scale && round(count1/y_scale) > count1/y_scale)
                    y = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                % before round value was closer to the left and lower neighbour
                elseif(round(count2/x_scale) > count2/x_scale && round(count1/y_scale) < count1/y_scale)
                    y = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) - 1,1),N0); % horizontal
                % before round value was closer to the right and lower neighbour
                else
                    y = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                end

                if(~isnan(B(y,x)))
                    Bmat(count1+1,count2+1) = B(y,x);
                else
                    Bmat(count1+1,count2+1) = 0;
                end
            end
        end

        % previous G1 positions
        for count1 = 0:2:M-1
            for count2 = 1:2:N-2
                % before round value was closer to the upper neighbour
                if(round(count1/y_scale) > count1/y_scale)
                    y = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale),1),N0); % horizontal
                else
                    y = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale),1),N0); % horizontal
                end

                if(~isnan(B(y,x)))
                    Bmat(count1+1,count2+1) = B(y,x);
                else
                    Bmat(count1+1,count2+1) = 0;
                end
            end
        end

        % previous G2 positions
        for count1 = 1:2:M-2
            for count2 = 0:2:N-1
                % before round value was closer to the left neighbour
                if(round(count2./x_scale) > count2/x_scale)
                    y = min(max(1+round(count1./y_scale),1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) - 1,1),N0); % horizontal
                else
                    y = min(max(1+round(count1./y_scale),1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                end

                if(~isnan(B(y,x)))
                    Bmat(count1+1,count2+1) = B(y,x);
                else
                    Bmat(count1+1,count2+1) = 0;
                end
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
                        Rmat(count1+1, count2+1) = mean(nonzeros(neighbors));
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
                        Gmat(count1+1, count2+1) = mean(nonzeros(neighbors));
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
                        Bmat(count1+1, count2+1) = mean(nonzeros(neighbors));
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

        % limit conditions
        Rmat(M,:) = Rmat(M-1,:);
        Rmat(:,N) = Rmat(:,N-1);

        Gmat(:,1) = Gmat(:,2);
        Gmat(M,:) = Gmat(M-1,:);
        Gmat(:,N) = Gmat(:,N-1);

        Bmat(M,:) = Bmat(M-1,:);
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
    
        % Split raw img to R, G, B channels
        R = NaN(M0,N0);
        G = NaN(M0,N0);
        B = NaN(M0,N0);

        R(Rx,Ry) = rawim(Rx,Ry);
        G(Gx1,Gy1) = rawim(Gx1,Gy1);
        G(Gx2,Gy2) = rawim(Gx2,Gy2);
        B(Bx,By) = rawim(Bx,By);

        % Create matrices to hold the interpolated R, G, B values
        Rmat = zeros(M,N);
        Gmat = zeros(M,N);
        Bmat = zeros(M,N);

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

        % Fill the rest of the R mat (previous G, B pixel spots) using nearest neighbor interpolation
        % previous G1 positions
        for count1 = 0:2:M-1
            for count2 = 0:2:N-1
                % before round value was closer to the left neighbour
                if(round(count2./x_scale) > count2/x_scale)
                    y = min(max(1+round(count1./y_scale),1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) - 1,1),N0); % horizontal
                else
                    y = min(max(1+round(count1./y_scale),1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                end

                if(~isnan(R(y,x)))
                    Rmat(count1+1,count2+1) = R(y,x);
                else
                    Rmat(count1+1,count2+1) = 0;
                end
            end
        end

        % previous G2 positions
        for count1 = 1:2:M-2
            for count2 = 1:2:N-2
                % before round value was closer to the upper neighbour
                if(round(count1/y_scale) > count1/y_scale)
                    y = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale),1),N0); % horizontal
                else
                    y = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale),1),N0); % horizontal
                end

                if(~isnan(R(y,x)))
                    Rmat(count1+1,count2+1) = R(y,x);
                else
                    Rmat(count1+1,count2+1) = 0;
                end
            end
        end

        % previous B positions
        for count1 = 1:2:M-2
            for count2 = 0:2:N-1
                % before round value was closer to the left and upper neighbour
                if(round(count2/x_scale) > count2/x_scale && round(count1/y_scale) > count1/y_scale)
                    y = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) -1,1),N0); % horizontal
                % before round value was closer to the right and upper neighbour
                elseif(round(count2/x_scale) < count2/x_scale && round(count1/y_scale) > count1/y_scale)
                    y = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                % before round value was closer to the left and lower neighbour
                elseif(round(count2/x_scale) > count2/x_scale && round(count1/y_scale) < count1/y_scale)
                    y = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) - 1,1),N0); % horizontal
                % before round value was closer to the right and lower neighbour
                else
                    y = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                end

                if(~isnan(R(y,x)))
                    Rmat(count1+1,count2+1) = R(y,x);
                else
                    Rmat(count1+1,count2+1) = 0;
                end
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
        for count1 = 0:2:M-1
            for count2 = 1:2:N-2
                % distances from borders: left, right, upper, lower 
                dist = [abs(round(count2/x_scale)-count2/x_scale), 1 - abs(round(count2/x_scale)-count2/x_scale), abs(round(count1/y_scale)-count1/y_scale), 1 - abs(round(count1/y_scale)-count1/y_scale)];
                [~,I] = min(dist);
                if(I == 1) % left is the nearest
                    y = min(max(1+round(count1./y_scale),1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) - 1,1),N0); % horizontal
                elseif(I == 2) % right
                    y = min(max(1+round(count1./y_scale),1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                elseif(I == 3) % upper
                    y = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale),1),N0); % horizontal
                else % lower
                    y = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale),1),N0); % horizontal
                end

                if(~isnan(G(y,x)))
                    Gmat(count1+1,count2+1) = G(y,x);
                else
                    Gmat(count1+1,count2+1) = 0;
                end
            end
        end

        % previous B positions
        for count1 = 1:2:M-2
            for count2 = 0:2:N-1
                % distances from borders: left, right, upper, lower 
                dist = [abs(round(count2/x_scale)-count2/x_scale), 1 - abs(round(count2/x_scale)-count2/x_scale), abs(round(count1/y_scale)-count1/y_scale), 1 - abs(round(count1/y_scale)-count1/y_scale)];
                [~,I] = min(dist);
                if(I == 1) % left is the nearest
                    y = min(max(1+round(count1./y_scale),1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) - 1,1),N0); % horizontal
                elseif(I == 2) % right
                    y = min(max(1+round(count1./y_scale),1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                elseif(I == 3) % upper
                    y = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale),1),N0); % horizontal
                else % lower
                    y = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale),1),N0); % horizontal
                end

                % if might not be necessary (2)
                if(~isnan(G(y,x)))
                    Gmat(count1+1,count2+1) = G(y,x);
                else
                    Gmat(count1+1,count2+1) = 0;
                end
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
        for count1 = 0:2:M-1
            for count2 = 1:2:N-2
                % before round value was closer to the left and upper neighbour
                if(round(count2/x_scale) > count2/x_scale && round(count1/y_scale) > count1/y_scale)
                    y = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) -1,1),N0); % horizontal
                % before round value was closer to the right and upper neighbour
                elseif(round(count2/x_scale) < count2/x_scale && round(count1/y_scale) > count1/y_scale)
                    y = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                % before round value was closer to the left and lower neighbour
                elseif(round(count2/x_scale) > count2/x_scale && round(count1/y_scale) < count1/y_scale)
                    y = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) - 1,1),N0); % horizontal
                % before round value was closer to the right and lower neighbour
                else
                    y = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                end

                if(~isnan(B(y,x)))
                    Bmat(count1+1,count2+1) = B(y,x);
                else
                    Bmat(count1+1,count2+1) = 0;
                end
            end
        end

        % previous G1 positions
        for count1 = 0:2:M-1
            for count2 = 0:2:N-1
                % before round value was closer to the upper neighbour
                if(round(count1/y_scale) > count1/y_scale)
                    y = min(max(1+round(count1./y_scale) - 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale),1),N0); % horizontal
                else
                    y = min(max(1+round(count1./y_scale) + 1,1),M0); % vertical
                    x = min(max(1+round(count2./x_scale),1),N0); % horizontal
                end

                if(~isnan(B(y,x)))
                    Bmat(count1+1,count2+1) = B(y,x);
                else
                    Bmat(count1+1,count2+1) = 0;
                end
            end
        end

        % previous G2 positions
        for count1 = 1:2:M-2
            for count2 = 1:2:N-2
                % before round value was closer to the left neighbour
                if(round(count2./x_scale) > count2/x_scale)
                    y = min(max(1+round(count1./y_scale),1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) - 1,1),N0); % horizontal
                else
                    y = min(max(1+round(count1./y_scale),1),M0); % vertical
                    x = min(max(1+round(count2./x_scale) + 1,1),N0); % horizontal
                end

                if(~isnan(B(y,x)))
                    Bmat(count1+1,count2+1) = B(y,x);
                else
                    Bmat(count1+1,count2+1) = 0;
                end
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
