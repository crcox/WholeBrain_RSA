function [X,info] = Adlas2(A,B,lambda, lambda1,options)
% Adlas  Sorted L1 parameter estimation solver
%
% [x,info] = ADLAS(A,b,lambda,options) solves the Slope problem
%
%     Minimize 1/2*||AX-B||_F^2 + sum_i OWL_lambda(x_i) + lambda1 ||X||_1,2
%
% where x_i denotes the i-th column in X. The entries in
% lambda must be nonnegative and in non-increasing order. 
%
% The options parameter is a structure with the following optional fields
% with [default value]: 
%
%   .iterations   Maximum number of iterations            [10,000]
%   .verbosity    0 = nothing, 1 = major, 2 = every             [1]
%   .fid        File identifier for output            [1 = stdout]
%   .optimIter    Iterations between optimality-condition checks    [1]
%   .gradIter     Iterations between full gradient computations    [20]
%   .tolInfeas    Maximum allowed dual infeasibility          [1e-6]
%   .tolRelGap    Stopping criterion for relative primal-dual gap [1e-6]
%   .xInit        Initial value of x                  [zeros(n,1)]
%
% The info output structure contains the following fields
%
%   .runtime      Runtime
%   .Aprods     Number of products with A
%   .ATprods      Number of products with A^T
%   .objPrimal    Primal objective
%   .objDual      Dual objective (possibly for infeasible dual point)
%   .infeas     Dual infeasibility
%   .status     Status: 1 = optimal, 2 = iterations
%

% Copyright 2013, M. Bogdan, E. van den Berg, W. Su, and E.J. Candes

% This file is part of SLOPE Toolbox version 1.0.
%
%   The SLOPE Toolbox is free software: you can redistribute it
%   and/or  modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of
%   the License, or (at your option) any later version.
%
%   The SLOPE Toolbox is distributed in the hope that it will
%   be useful, but WITHOUT ANY WARRANTY; without even the implied
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%   See the GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with the SLOPE Toolbox. If not, see
%   <http://www.gnu.org/licenses/>.

  % -------------------------------------------------------------
  % Start timer
  % -------------------------------------------------------------
  t0 = tic();


  % -------------------------------------------------------------
  % Parse parameters
  % -------------------------------------------------------------
  if (nargin <  4), options = struct(); end;

  miniterations = getDefaultField(options,'miniterations',1000);
  iterations = getDefaultField(options,'iterations',10000);
  verbosity  = getDefaultField(options,'verbosity',0);
  fid     = getDefaultField(options,'fid',1);
  optimIter  = getDefaultField(options,'optimIter',1);
  gradIter   = getDefaultField(options,'gradIter',20);
  tolInfeas  = getDefaultField(options,'tolInfeas',1e-7);
  tolRelGap  = getDefaultField(options,'tolRelGap',1e-7);
  xInit     = getDefaultField(options,'xInit',[]);

  % Ensure that lambda is non-increasing
  if ((length(lambda) > 1) && any(lambda(2:end) > lambda(1:end-1)))
    error('Lambda must be non-increasing.');
  end
  if (lambda(end) < 0)
    error('Lambda must be nonnegative');
  elseif (lambda(1) == 0)
    error('Lambda must have at least one nonnegative entry.');
  end


  % -------------------------------------------------------------
  % Initialize
  % -------------------------------------------------------------

  % Get problem dimension
  n = size(A,2);
  r = size(B,2);
  % Get initial lower bound on the Lipschitz constant
  s = RandStream('mt19937ar','Seed',0);
  X = randn(s,n,r); X = X / norm(X,'fro');
  X = A'*(A*X);
  L = norm(X,'fro');
  L = 1;

  % Constants for exit status
  STATUS_RUNNING    = 0;
  STATUS_OPTIMAL    = 1;
  STATUS_ITERATIONS = 2;
  STATUS_ALLZERO    = 3;
  STATUS_MSG = {'Optimal','Iteration limit reached','All weights set to zero'};

  % Initialize parameters and iterates
  if (isempty(xInit)), xInit = zeros(n,r); end;

  t     = 1;
  eta   = 2;
  lambda  = lambda(:);
  X     = xInit;
  Y     = X;
  Ax      = A*X;
  fPrev   = Inf;
  iter    = 0;
  status  = STATUS_RUNNING;
  Aprods  = 2;
  ATprods = 1;
  U = zeros(size(X));
  Z = X;

  % Deal with Lasso case
  modeLasso = (numel(lambda) == 1);
  if (modeLasso)
    proxFunction = @(v1,v2) proxL1L2(v1,v2);
  else
    proxFunction = @(v1,v2) proxGrowl2(v1,v2);
  end

  if (verbosity > 0)
    fprintf(fid,'%5s  %9s %9s  %9s  %9s\n','Iter','||r||_F','Gap','Infeas.','Rel. gap');
  end


  % -------------------------------------------------------------
  % Main loop
  % -------------------------------------------------------------
  while (true)

    % Compute the gradient at f(y)
    if (mod(iter,gradIter) == 0) % Includes first iterations
      r = A*Y - B;
      g = A'*(A*Y-B);
      f = trace(r'*r) / 2;
    else
      r = (Ax + ((tPrev - 1) / t) * (Ax - AxPrev)) - B;
      g = A'*(A*Y-B);
      f = trace(r'*r) / 2;
    end

    % Increment iteration count
    iter = iter + 1;

    % Check optimality conditions
    if ((mod(iter,optimIter) == 0))
      if iter>=miniterations
        if (abs( fPrev - f ) < tolRelGap*max(fPrev,1))
          status = STATUS_OPTIMAL;
          break;
        end
      end
    else
      str = '';
    end

    if (verbosity > 0)
      if ((verbosity == 2) || ...
        ((verbosity == 1) && (mod(iter,optimIter) == 0)))
      fprintf(fid,'%5d  %9.2e%s\n', iter,f,str);
      end
    end

    % Stopping criteria
    if (status == 0)
      if (iter >= iterations)
        status = STATUS_ITERATIONS;
      end
    end

    if (status ~= 0)
      if (verbosity > 0)
        fprintf(fid,'Exiting with status %d -- %s\n', status, STATUS_MSG{status});
      end
      break;
    end

    % Keep copies of previous values
    AxPrev = Ax;
    xPrev  = Z;%X;
    fPrev  = f;
    tPrev  = t;

    % Lipschitz search
    while (L < inf)
      % Compute prox mapping

      X = proxFunction( (Z - U) - (1/L)*g, lambda/L);

      d = X - Y;
      %disp(L)
      Ax = A*X;%A1*vec(X);
      r = Ax-B;
      f = trace(r'*r)/2;
      q = fPrev + trace(d'*g) + (L/2)*trace(d'*d);

      Aprods = Aprods + 1;

      if (q >= f*(1-1e-12))
        break;
      else
        L = L * eta;
      end
    end
    % Update
    t = (1 + sqrt(1 + 4*t^2)) / 2;
    Y = X + ((tPrev - 1) / t) * (X - xPrev);

    Z = proxL1L2(Y + U, lambda1/L);
    U = U + (Y - Z);

    % Check if all weights set to zero
    if all(Y(:)==0) && iter > 100
      status = STATUS_ALLZERO;
      break
    end
  end

  % Set solution
  X = Y;
  %eps = max(sqrt(sum(X.^2,2)))*0.05;
  %X(sqrt(sum(X.^2,2))<0.001,:) = 0;
  % Information structure
  info = struct();
  if (nargout > 1)
    info.runtime   = toc(t0);
    info.Aprods    = Aprods + ceil(iter / gradIter);
    info.ATprods   = ATprods + iter;
    info.objPrimal = nan(1);
    info.objDual   = nan(1);
    info.infeas    = nan(1);
    info.status    = status;
    info.message   = STATUS_MSG{status};
    info.iter      = iter;
    info.L         = L;
  end
end % Function Adlas


% ------------------------------------------------------------------------
function opt = getDefaultField(data,field,default)
% ------------------------------------------------------------------------
  if isfield(data,field)
    opt = data.(field);
  else
    opt = default;
  end
end


function x = proxL1L2(Y,lambda)
% ------------------------------------------------------------------------
  % Normalization
    %y    = sqrt(sum(Y.^2,2));
  %disp('here')
   tmp = Y;
  r = size(Y,2);
  xtmp = tmp./(repmat(sqrt(sum(tmp.^2,2))+realmin,1,r));
   x = xtmp.*repmat(max(sqrt(sum(tmp.^2,2))-lambda,0),1,r);
  %disp(nnz(any(x,2)))
  %x    = Y .* ((max(y - lambda,0)./(y + realmin))*ones(1,size(Y,2)));
end

function X = proxGrowl2(Y, lambda)
   %run OWL prox on each column separately
   lambda = lambda(:);
   X = zeros(size(Y));

   r = size(Y,2);
   for ii = 1:r
     X(:,ii) = proxSortedL1(Y(:,ii),lambda);
   end
end
