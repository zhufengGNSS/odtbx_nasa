classdef estimator < handle
    % ESTIMATOR Estimator base class
    %   This class is used to define the properties common to all ODTBX
    %   estimators. It will promote uniformity in implementations and
    %   coding.
    
    properties
        t
        Xhat
        Phat
        e
        y
        Pa
        Pv
        Pw
        Phata
        Phatv
        Phatw
        Sigma_a % estbat
        Sig_sa % estseq
        eflag % estseq
        Pdy
        Pdyt
        Pm % estseq
        Phatm % estseq
        restartRecord % estseq
    end
    
    methods
        function obj = estimator()
            % Preallocate these variables are cells (makes it faster?)
            obj.t = {}; 
            obj.Xhat = {};
            obj.Phat = {};
            obj.e = {};
            obj.y = {};
            obj.Pdy = {};
            obj.eflag = {};
        end
        
        function varargout = run_estimator(obj,varargin)
            % This is a dummy function. Subclasses of estimator should
            % overwrite this function with the function that will run the
            % estimator.
        end
    end % Methods
    
    methods(Static)
        function x = robustls(A,b)
            % More robust least-squares solution to Ax = b.  This is based on the help
            % for the QR function, which shows how "the least squares approximate
            % solution to A*x = b can be found with the Q-less qr decomposition and one
            % step of iterative refinement."
            R = triu(qr(A));
            x = R\(R'\(A'*b));
            r = b - A*x;
            e = R\(R'\(A'*r));
            x = x + e;
        end

        % Self-test user functions
        function [Xdot,A,Q] = rwdyn(t,X,q) % Test 1 dynfun
            el = length(t);
            A = zeros(1,1,el);
            Xdot = A*X;
            Q = q*ones(1,1,el);
        end

        function [Y,H,R] = rwdat(t,X,r) % Test 1 datfun
            Y = X;
            H = ones(1,1,length(t));
            R = r*ones(1,1,length(t));
        end

        function [Xdot,A,Q] = irwbdyn(t,X,q) % Test 2 dynfun.tru
            el = length(t);
            A([1:3 7:8],[4:6 8]) = eye(5,4);
            Xdot = A*X;
            A = repmat(A,[1 1 el]);
            Q([4:6 8],[4:6 8]) = q*eye(4);
            Q = repmat(Q,[1 1 el]);
        end
        
        function [Xdot,A,Q] = irwdyn(t,X,q) % Test 2 dynfun.est
            el = length(t);
            A(:,4:6) = eye(6,3);
            Xdot = A*X;
            A = repmat(A,[1 1 el]);
            Q(4:6,4:6) = q*eye(3);
            Q = repmat(Q,[1 1 el]);
        end
        
        function [Xdot,A,Q] = sogmbdyn(t,X,c) % Test 3 dynfun.tru
            el = length(t);
            A(:,2:3) = eye(3,2);
            A(3,2:3) = [-c.w_n^2, -2*c.zeta*c.w_n];
            Xdot = A*X;
            A = repmat(A,[1 1 el]);
            Q(3,3) = c.q;
            Q = repmat(Q,[1 1 el]);
        end
        
        function [Y,H,R] = sogmbdat(t,X,r) % Test 3 datfun.tru
            Y = X(1,:) + X(2,:);
            H = [1 1 0];
            H = repmat(H,[1,1,length(t)]);
            R = r*ones(1,1,length(t));
        end
        
        function [Y,H,R] = irwbdat(t,X,r) % Test 2 datfun.tru
            el = length(t);
            H = [eye(3) zeros(3,3) -ones(3,1) zeros(3,1)];
            Y = H*X;
            H = repmat(H,[1 1 el]);
            R = r*repmat(eye(3),[1 1 el]);
        end
        
        function [Y,H,R] = irwdat(t,X,r) % Test 2 datfun.est
            el = length(t);
            H = [eye(3) zeros(3,3)];
            Y = H*X;
            H = repmat(H,[1 1 el]);
            R = r*repmat(eye(3),[1 1 el]);
        end
        
        function y = refine(u,refine)
            y = [reshape([u(1:end-1);repmat(u(1:end-1),refine,1)+...
            cumsum(repmat(diff(u)/(refine+1),refine,1),1)],[],1);u(end)]';
        end
    end % Methods (Static)
end % Class

