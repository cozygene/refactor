
refactor <- function(data_file, k, t=500, numcomp=NULL, ranked_filename='refactor.out.rankedlist.txt', components_filename='refactor.out.components.txt') {

    print('Starting ReFACTor v1.0...');

    print('Reading input files...');

    O = as.matrix(read.table(data_file))
    sample_id <- O[1, -1] # extract samples ID
    O <- O[-1,] # remove sample ID from matrix
    cpgnames <- O[, 1] ## set rownames
    O <- O[, -1] 
    O = matrix(as.numeric(O),nrow=nrow(O),ncol=ncol(O))

    if (is.null(numcomp) || is.na(numcomp)) 
    {
        numcomp = k
    }
    
    print('Running a standard PCA...')
    pcs = prcomp(scale(t(O)));

    coeff = pcs$rotation
    score = pcs$x

    print('Compute a low rank approximation of input data and rank sites...')
    x = score[,1:k]%*%t(coeff[,1:k]);
    An = scale(t(O),center=T,scale=F)
    Bn = scale(x,center=T,scale=F)
    An = t(t(An)*(1/sqrt(apply(An^2,2,sum))))
    Bn = t(t(Bn)*(1/sqrt(apply(Bn^2,2,sum))))


    # Find the distance of each site from its low rank approximation.
    distances = apply((An-Bn)^2,2,sum)^0.5 ;
    dsort = sort(distances,index.return=T);
    ranked_list = dsort$ix

    print('Compute ReFACTor components...')
    sites = ranked_list[1:t];
    pcs = prcomp(scale(t(O[sites,])));
    first_score <- score[,1:k];
    score = pcs$x

    print('Saving a ranked list of the data features...');
    write(t(cbind(ranked_list,cpgnames[ranked_list])),file=ranked_filename,ncol=2)

    print('Saving the ReFACTor components...');
    write(t(score[,1:numcomp]), file=components_filename, ncol=numcomp)
    
    print('ReFACTor is Done');
    result <- list(refactor_components=score[,1:numcomp], ranked_list=ranked_list, standard_pca=first_score) 
    return(result)

}
