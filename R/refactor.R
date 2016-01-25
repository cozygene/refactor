
refactor <- function(data_file, k, covarfile = NULL, t = 500, numcomp = NULL, stdth = 0.02, out = "refactor") {

    ranked_filename = paste(out, ".out.rankedlist.txt", sep="")
    components_filename = paste(out, ".out.components.txt", sep="")

    print('Starting ReFACTor v1.0...');

    print('Reading input files...');

    O = as.matrix(read.table(data_file))
    sample_id <- O[1, -1] # extract samples ID
    O <- O[-1,] # remove sample ID from matrix
    cpgnames <- O[, 1] ## set rownames
    O <- O[, -1] 
    O = matrix(as.numeric(O),nrow=nrow(O),ncol=ncol(O))

    print(paste("Excluding sites with low variance (std < ", stdth, ")..."), sep="")
    sds = apply(t(O), 2, sd)
    m_before = length(sds)
    include = which(sds >= stdth)
    O = O[include,]
    cpgnames = cpgnames[include]
    print(paste((m_before - length(which(sds >= stdth))), " sites were excluded due to low variance...", sep=""))

    if (is.null(numcomp) || is.na(numcomp)) 
    {
        numcomp = k
    }

    # Adjust the data for the covariates
    if (!is.null(covarfile))
    {
        covs = as.matrix(read.table(covarfile))
        sample_id2 <- covs[, 1]
        if (!all(sample_id == sample_id2)){
            print("ERROR: The order of the samples in the covariates file must be the same as the order in the data file")
            quit()
        }
        covs <- covs[,-1]
        if (length(covs) > dim(O)[2])
        {
            covs = matrix(as.numeric(covs),nrow=nrow(covs),ncol=ncol(covs))
        }else{
            covs = as.numeric(covs)
        }    
        O_adj = O
        for (site in 1:nrow(O))
        {
            model <- lm(O[site,] ~  covs)
            O_adj[site,] = residuals(model)
        }
        O = O_adj
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
    write(t(cpgnames[ranked_list]),file=ranked_filename,ncol=1)
    #write(t(cbind(ranked_list,cpgnames[ranked_list])),file=ranked_filename,ncol=2)

    print('Saving the ReFACTor components...');
    write(t(score[,1:numcomp]), file=components_filename, ncol=numcomp)
    
    print('ReFACTor is Done');
    result <- list(refactor_components=score[,1:numcomp], ranked_list=ranked_list, standard_pca=first_score) 
    return(result)

}
