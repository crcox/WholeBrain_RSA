function predictedSimilarity(result,varargin)
    p = inputParser();
    addParameter(p,'filter','all',@ischar);
    parse(p, varargin{:});
    
    filter = p.Results.filter;
    cvfilter = result.cvfilter(~result.finalfilter);
    
    switch filter
        case 'all'
            figure
            subplot(1,2,1)
            imagesc(result.S);
            title('True Similarity')
            subplot(1,2,2)
            imagesc(result.Sz);
            title('Predicted Similarity');
        case 'test'
            figure
            subplot(1,2,1)
            imagesc(result.S(cvfilter,cvfilter));
            title('True Similarity')
            subplot(1,2,2)
            imagesc(result.Sz(cvfilter,cvfilter));
            title('Predicted Similarity');
        case 'train'
            figure
            subplot(1,2,1)
            imagesc(result.S(~cvfilter,~cvfilter));
            title('True Similarity')
            subplot(1,2,2)
            imagesc(result.Sz(~cvfilter,~cvfilter));
            title('Predicted Similarity');
    end
end
