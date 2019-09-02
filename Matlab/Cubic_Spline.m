classdef Cubic_Spline
    % Cubic Spline class allows for interpolation between data points
    % 
    % @author: Matt Marti
    % @date: 2019-07-20
    
    properties
        akmat
        bkmat
        ckmat
        dkmat
        xkvec
        extrap_flag = 1;
    end
    
    methods
        function this = Cubic_Spline(xkvec, fkvec, fslope)
            % Cubic spline interpolation function
            % Solves for the parameters that describe a cubic spline. This 
            % function can solve the parameters of a cubic spline given 
            % data from a vectorized system. For example, the three space 
            % position of a planet in orbit can be interpolated.
            % 
            % @arg
            % xkvec         - 1 x n double matrix
            %                 Independent variable data points
            % fkvec         - m x n double matrix
            %                 Dependent variable data points
            % fslope        - m x 2 double matrix (optional)
            %                 Function slope at boundary points
            % 
            % @return
            % splineDataMat - m x n x 5 double matrix
            %                 Spline coefficient data matrix. Organized by 
            %                 input data dimension, known value points, and
            %                 coefficient.
            % 
            % @author: Matt Marti
            % @date: 2019-07-21
            
            
            % Size
            nx = length(xkvec);
            m = size(fkvec, 1);

            % Check input array sizes
            assert(size(xkvec,1) == 1 || size(xkvec,2) == 1, ...
                'Argument ''xkvec'' is not a vector');
            if size(xkvec,1) ~= 1
                xkvec = xkvec';
            end
            this.xkvec = xkvec;
            assert(size(fkvec,2) == nx, ...
                'Argument ''fkvec'' does not match length of ''xkvec''');
            if nargin > 2
                assert(size(fslope,1) == m && size(fslope,2) == 2, ...
                    'Argument ''fslope'' is of incorrect size: is %d by %d', ...
                    size(fslope,1), size(fslope,2));
                cbcflag = 1;
            else
                cbcflag = 0;
            end

            % Build tri-diagonal system of equations
            hkvec = xkvec(2:nx) - xkvec(1:nx-1);
            n = nx - 1;
            H = eye(nx);
            for k = 2:n
                H(k, k-1) = hkvec(k-1);
                H(k, k) = 2 * ( hkvec(k-1) + hkvec(k) );
                H(k, k+1) = hkvec(k);
            end
            if cbcflag
                H(1, 1) = 2*hkvec(1);
                H(1, 2) = hkvec(1);
                H(n+1, n) = hkvec(n);
                H(n+1, n+1) = 2*hkvec(n);
            end

            % Generate right hand side
            this.akmat = fkvec';
            xstar = zeros(nx,m);
            for k = 2:n
                xstar(k,:) = 3*( ( this.akmat(k+1,:) - this.akmat(k,:) ) / hkvec(k) ...
                              - ( this.akmat(k,:) - this.akmat(k-1,:) ) / hkvec(k-1) )';
            end
            if cbcflag
                xstar(1,:) = 3*( (this.akmat(2,:) - this.akmat(1,:) )/hkvec(1) - fslope(:,1)' );
                xstar(nx,:) = 3*( fslope(:,2)' - (this.akmat(nx,:) - this.akmat(nx-1,:) )/hkvec(nx-1) );
            end

            % Solve tri-diagonal system of equations
            this.ckmat = H \ xstar;

            % Compute bkmat and dkmat
            this.bkmat = zeros(nx,m);
            this.dkmat = zeros(nx,m);
            for k = 1:n
                this.bkmat(k,:) = ( this.akmat(k+1,:) - this.akmat(k,:) ) / hkvec(k) ...
                    - hkvec(k) * ( 2*this.ckmat(k,:) + this.ckmat(k+1,:) ) / 3;
                this.dkmat(k,:) = ( this.ckmat(k+1,:) - this.ckmat(k,:) ) / ( 3*hkvec(k) );
            end
            if cbcflag
                this.bkmat(nx,:) = fslope(:,2)';
            end
        end
        
        function [finter, dfinter] = interp(this, xinter, dflag)
            % Cubic spline interpolation function
            % Interpolates function values at specified points using data 
            % from a solved cubic spline.
            % 
            % @arg
            % splineDataMat - m x n x 5 double matrix
            %                 Spline coefficient data matrix. Organized by 
            %                 input data dimension, known value points, and
            %                 coefficient.
            % xinter        - 1 x n double matrix
            %                 Interpolation points
            % dflag         - bool (optional)
            %                 Optional flag to make the function return the
            %                 derivative.
            %                 False by default.
            % 
            % @return
            % finter        - n x 1 double matrix
            %                 Interpolated function value
            % dfinter       - n x 1 double matrix
            %                 Interpolated function derivative value
            % 
            % @author: Matt Marti
            % @date: 2019-07-21
            
            % Check input
            if nargin < 3
                dflag = 0;
            end
            assert(size(xinter,1) == 1 || size(xinter,2) == 1, ...
                'Argument ''xinter'' is not a vector');

            % Size
            nx = size(this.akmat, 1);
            m = size(this.akmat, 2);
            
            % Preallocate output
            finter = zeros(m,length(xinter));
            if dflag
                dfinter = zeros(m,length(xinter));
            else
                dfinter = [];
            end
            
            % Interpolate function
            for j = 1:length(xinter)

                % Check that interpolated value is within function range
                if ~this.extrap_flag
                    assert(xinter(j) >= this.xkvec(1) && xinter(j) <= this.xkvec(nx),...
                        'Interpolation value not within bounds');
                end

                % Find x value just below xinter(i)
                k = 2;
                while k <= nx
                    if xinter(j) < this.xkvec(k), k = k - 1; break; end
                    k = k + 1;
                end
                if k > nx % Point is on upper boundary
                    finter(:,j) = this.akmat(nx,:);
                    if dflag
                        dfinter(:,j) = this.bkmat(nx,:);
                    end
                    continue;
                end

                % Spline interpolation
                hi = xinter(j) - this.xkvec(k);
                finter(:,j) = this.akmat(k,:) + this.bkmat(k,:)*hi + this.ckmat(k,:)*hi^2 + this.dkmat(k,:)*hi^3;
                if dflag
                    dfinter(:,j) = this.bkmat(k,:) + 2*this.ckmat(k,:)*hi + 3*this.dkmat(k,:)*hi^2;
                end
            end
        end
    end
end