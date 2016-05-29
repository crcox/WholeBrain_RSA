function [Aredux, reduxFilter]=removeOutliers(A)
% This function will z-score the row and column means and remove those
% which have a z-score greater than 5.  The process it performed on
% columns, and then rows.  Each time outliers are removed, the reduced data
% is checked again for outliers, until no outliers remain.

    reduxFilter.words=true(1,size(A,1));
    reduxFilter.voxels=true(1,size(A,2));
    index.row = 1:size(A,1);
    index.col = 1:size(A,2);
    Aredux=A;

    % For columns...
    outVector = abs(zscore(mean(Aredux,1))) > 5;
    outCount = sum(outVector);
    outliers = [];
    while outCount > 0
        outliers = [outliers index.col(outVector)]; %#ok<AGROW>
        index.col = index.col(~outVector);
        Aredux = Aredux(:,~outVector);
        outVector = abs(zscore(mean(Aredux,1))) > 5;
        outCount = sum(outVector);
    end
    reduxFilter.voxels(outliers) = false;
    reduxFilter.voxels = reduxFilter.voxels & any(A,1);

    % For rows...
    outVector = abs(zscore(mean(Aredux,2))) > 5;
    outVector = outVector'; % Mean by row creates a vector n-by-1 instead of 1-by-n.
    outCount = sum(outVector);
    outliers = [];
    while outCount > 0
        outliers = [outliers index.row(outVector)]; %#ok<AGROW>
        index.row = index.row(~outVector);
        Aredux = Aredux(~outVector,:);
        outVector = abs(zscore(mean(Aredux,2))) > 5;
        outVector = outVector'; % Mean by row creates a vector n-by-1 instead of 1-by-n.
        outCount = sum(outVector);
    end
    reduxFilter.words(outliers) = false;
    reduxFilter.words = reduxFilter.words & any(A,2)';
end
