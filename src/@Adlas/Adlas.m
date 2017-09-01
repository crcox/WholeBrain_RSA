classdef Adlas
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        A
        B
        lambda
        X % model weights
        xInit
        max_iter = 100000
        fid = 1
        optimIter = 1
        gradIter = 20
        tolInfeas = 1e-6
        tolRelGap = 1e-8
        n % size(A,2)
        r % size(B,2)
        s % RandStream('mt19937ar','Seed',0)
        L % Lipschitz constant
        t = 1
        eta = 2
        iter = 0;
        status
        message
        verbosity = 0
        objPrimal
        objDual
        infeas
        Aprods
    end

    methods
        function obj = Adlas(A,B,lambda,opts)
            if (nargin <  4), opts = struct(); end;
            obj.A = A;
            obj.B = B;
            % Ensure that lambda is non-increasing
            if ((length(lambda) > 1) && any(lambda(2:end) > lambda(1:end-1)))
                error('Lambda must be non-increasing.');
            end
            if (lambda(end) < 0)
                error('Lambda must be nonnegative');
            elseif (lambda(1) == 0)
                error('Lambda must have at least one nonnegative entry.');
            end
            obj.lambda = lambda;
            for i = 1:numel(fn)
                obj.(fn{i}) = opts.(fn{i});
            end
            obj.n = size(A,2);
            obj.r = size(B,2);
            obj.s = RandStream('mt19937ar','Seed',0);
            obj.L = 1;
            if ~isfield(opts, 'xInit') || (isempty(opts.xInit))
                xInit = zeros(obj.n,obj.r);
            end
        end
        function obj = train(opts)
            fn = fieldnames(opts);
            for i = 1:numel(fn)
                obj.(fn{i}) = opts.(fn{i});
            end
            obj = Adlas1(obj);
        end
    end
end

function obj = Adlas1(obj)
    % Constants for exit status
    STATUS_RUNNING    = 0;
    STATUS_OPTIMAL    = 1;
    STATUS_ITERATIONS = 2;
    STATUS_ALLZERO    = 3;
    STATUS_MSG = {'Optimal','Iteration limit reached','All weights set to zero'};

    % Initialize parameters
    X       = obj.xInit;
    Y       = X;
    Ax      = obj.A * X;
    status  = STATUS_RUNNING;

    % Deal with Lasso case
    modeLasso = (numel(lambda) == 1);
    if (modeLasso)
        proxFunction = @(v1,v2) proxL1L2(v1,v2);
    else
        proxFunction = @(v1,v2) proxSortedL1L2(v1,v2);
    end

    if (verbosity > 0)
        fprintf(fid,'%5s  %9s %9s  %9s  %9s\n','Iter','||r||_F','Gap','Infeas.','Rel. gap');
    end

    % -------------------------------------------------------------
    % Main loop
    % -------------------------------------------------------------
    while (true)

        % Compute the gradient at f(y)
        if (mod(obj.iter,gradIter) == 0) % Includes first iterations
            r = A*Y - B;
            g = A'*(A*Y-B);
            f = trace(r'*r) / 2;
        else
            r = (Ax + ((tPrev - 1) / obj.t) * (Ax - AxPrev)) - B;
            g = A'*(A*Y-B);
            f = trace(r'*r) / 2;
        end

        % Increment iteration count
        obj.iter = obj.iter + 1;

        % Check optimality conditions
        if ((mod(obj.iter,optimIter) == 0))
            % Compute 'dual', check infeasibility and gap
            if (modeLasso)
                gs = sqrt(sum(g.^2,2));
                ys = sqrt(sum(Y.^2,2));

                infeas = max(norm(gs,inf)-lambda,0);

                objPrimal = f + lambda*norm(ys,1);
                objDual   = -f - trace(r'*B);
            else
                gs     = sort(sqrt(sum(g.^2,2)),'descend');
                ys     = sort(sqrt(sum(Y.^2,2)),'descend');
                infeas = max(max(cumsum(gs-lambda)),0);

                % Compute primal and dual objective
                objPrimal =  f + lambda'*ys;
                objDual  = -f - trace(r'*B);
            end

            % Format string
            if (verbosity > 0)
                str = sprintf(' %9.2e  %9.2e  %9.2e',objPrimal - objDual, infeas/lambda(1), abs(objPrimal - objDual) / max(1,objPrimal));
            end

            % Check primal-dual gap
            if ((abs(objPrimal - objDual)/max(1,objPrimal) < tolRelGap)  && ...
                    (infeas < tolInfeas * lambda(1)) )
                status = STATUS_OPTIMAL;
            end

        else
            str = '';
        end

        if (verbosity > 0)
            if ((verbosity == 2) || ...
                    ((verbosity == 1) && (mod(obj.iter,optimIter) == 0)))
                fprintf(fid,'%5d  %9.2e%s\n', obj.iter,f,str);
            end
        end

        % Stopping criteria
        if (status == 0)
            if (obj.iter >= max_iter)
                status = STATUS_ITERATIONS;
            end
        end

        if (status ~= 0)
            if verbosity > 0
                fprintf(fid,'Exiting with status %d -- %s\n', status, STATUS_MSG{status});
            end
            break;
        end

        % Keep copies of previous values
        AxPrev = Ax;
        xPrev  = X;
        fPrev  = f;
        tPrev  = obj.t;

        % Lipschitz search
        while (obj.L < inf)
            % Compute prox mapping
            X = proxFunction(Y - (1/obj.L)*g, lambda/obj.L);
            d = X - Y;

            Ax = A*X;%A1*vec(X);
            r = Ax-B;
            f = trace(r'*r)/2;
            q = fPrev + trace(d'*g) + (obj.L/2)*trace(d'*d);

            obj.Aprods = obj.Aprods + 1;

            if (q >= f*(1-1e-12))
                break;
            else
                obj.L = obj.L * obj.eta;
            end
        end

        % Update
        obj.t = (1 + sqrt(1 + 4*obj.t^2)) / 2;
        Y = X + ((tPrev - 1) / obj.t) * (X - xPrev);

        % Check is all weights set to zero
        if all(Y(:)==0) && obj.iter > 100
            status = STATUS_ALLZERO;
            break
        end
    end

    % Set solution
    obj.X = Y;
    obj.objPrimal = objPrimal;
    obj.objDual   = objDual;
    obj.infeas    = infeas;
    obj.status    = status;
    obj.message   = STATUS_MSG{status};
    obj.Aprods   = obj.Aprods + ceil(obj.iter / gradIter);
end

function x = proxL1L2(Y,lambda)
    tmp = Y;
    r = size(Y,2);
    xtmp = tmp./(repmat(sqrt(sum(tmp.^2,2))+realmin,1,r));
    x = xtmp.*repmat(max(sqrt(sum(tmp.^2,2))-lambda,0),1,r);
end
