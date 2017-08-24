function [U, obj] = grouplasso_glmnet(X, Y, alpha, lambda, varargin)
    p = inputParser;
    %% Parse function inputs
    %                name         default     validation
    addRequired(p  , 'X'                                 );
    addRequired(p  , 'Y'                                 );
    addRequired(p  , 'alpha'                             );
    addRequired(p  , 'lambda'                            );
    addParameter(p , 'cvind'   , []                      );
    addParameter(p , 'bias'    , 0          , @isscalar  );
    addParameter(p , 'maxiter' , 1000       , @isscalar  );
    addParameter(p , 'tol'     , 1e-8       , @isscalar  );
    addParameter(p , 'U0'      , []                      );
    addParameter(p , 'PARALLEL', false                   );
    addParameter(p , 'verbose' , false                   );
    parse(p, X, Y, alpha, lambda, varargin{:});

    X         = p.Results.X;
    Y         = p.Results.Y;
    alpha     = p.Results.alpha;
    lambda    = p.Results.lambda;
    bias      = p.Results.bias;
    maxiter   = p.Results.maxiter;
    tol       = p.Results.tol;
    U0        = p.Results.U0;
    CVIND     = p.Results.cvind;
    PARALLEL   = p.Results.PARALLEL;
    verbose   = p.Results.verbose;

    if ~iscell(X)
        X = {X};
    end
    if ~iscell(Y)
        Y = {Y};
    end
    if ~iscell(CVIND)
        CVIND = {CVIND};
    end

    U = cell(numel(X),1);
    for iSubj = 1:numel(X)
        x = X{iSubj};
        y = Y{iSubj};

        cvind = CVIND{iSubj};
        cvset = unique(cvind);
        nfold = numel(cvset);

        cinds = unique(y);
        cinds = cinds(:)'; % force to row vec
        m = numel(cinds);

        performanceMetric = 'mse';
        modelType = 'mgaussian';

        if cvset(1) > 1
            cvind = cvind - (cvset(1)-1);
            cvset = cvset - (cvset(1)-1);
        end
        adj = find(diff(cvset)>1);
        if ~isempty(adj);
            for a = adj
                cvind(cvind>a) = cvind(cvind>a) - 1;
            end
        end
        if isnan(lambda)
            opts = glmnetSet(struct('intr',bias,'thresh', tol, 'weights', U0, 'alpha', alpha, 'lambda', []));
            objcv = cvglmnet(x,y,modelType,opts,performanceMetric,nfold,cvind',PARALLEL);
            lambda = objcv.lambda_min;
        end
        opts = glmnetSet(struct('intr',bias,'thresh', tol, 'weights', U0, 'alpha', alpha, 'lambda', lambda));
        obj(iSubj) = glmnet(x,y,modelType,opts);
        if iscell(obj(iSubj).beta)
            U{iSubj} = cell2mat(obj(iSubj).beta);
        else
            U{iSubj} = obj(iSubj).beta;
        end
    end
end
