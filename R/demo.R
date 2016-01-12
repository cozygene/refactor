
if(!exists("refactor", mode="function")) source("refactor.R")

associations_test <- function(O, y, model_append)
{
    observed_pvalues <- c()
    m <- ncol(O)
    for (site in 1:m)
    {
        model <- lm(y ~ O[,site] + model_append)
        
        pvalue <- coef(summary(model))[2,4]
        observed_pvalues[site] = as.numeric(pvalue)
    }

    x <- runif(m)
    x <- sort(-log10(x))
    y <- sort(-log10(observed_pvalues))
    result <- list(x=x,y=y) 
    return(result) 
}


args <- commandArgs(trailingOnly = TRUE)
data_file=args[1] # data file
y=args[2] # a phenotype file
K=as.numeric(args[3]) # an integer
R=args[4] # a matrix with k columns; these are the cell proportions)

## TODO check that arguments are not NULL
# if (data_file == NULL) {
#     print("Running on data file: demo/demo_datafile.txt")
#     data_file <- "demo/demo_datafile.txt"
# }
# if (y == NULL) {
#     print("Running on phenotype file: demo/demo_phenotype.txt")
#     y <- "demo/demo_phenotype.txt"
# }
# if (K == NULL) {
#     print("Running on K = 5")
#     K<-5
# }
# if (R == NULL) {
#     print("Running on cellproportions file: demo/demo_cellproportions.txt")
#     R <- "demo/demo_cellproportions.txt"
# }  

#TODO enter to funs?
    O = as.matrix(read.table(data_file))    
    sample_id_O <- O[1, -1] # extract samples ID
    O <- O[-1,] # remove sample ID from matrix
    cpgnames <- O[, 1] ## set rownames
    O <- O[, -1] 
    O = t(matrix(as.numeric(O),nrow=nrow(O),ncol=ncol(O)))

phenotype_matrix = as.matrix(read.table(y))
y <-  matrix(as.numeric(as.matrix(phenotype_matrix[, -1]) ))

sample_id_y <- phenotype_matrix[, 1]

# is there a better way? this takes time
for (i in 1:length(sample_id_y)){
    if (sample_id_O[i] != sample_id_y[i]) {
        print( "ERROR IN SAMPLE ID")
    }
}


png('plot.png')
par(mfrow=c(2,2))

#for printing the x=y trend
xasix<-0:10
yasix<-0:10


output <- refactor(data_file, K)

# exp1
print("START EXP1")
res <- associations_test(O,y, matrix(0, nrow = nrow(O), ncol = 1))
plot(res$x, res$y, main="original", xlab="uniform distribution", ylab="observed pvalues", pch='.', xlim=c(0,ceiling(max(res$x))),  ylim=c(0,ceiling(max(res$y))))
abline(lm(xasix~yasix),col=2,lty=1)


# exp2
print("START EXP2")
res <- associations_test(O, y, as.matrix(read.table(R)))
plot(res$x, res$y, main="R", xlab="uniform distribution", ylab="observed pvalues", pch='.', xlim=c(0,ceiling(max(res$x))),  ylim=c(0,ceiling(max(res$y))))
abline(lm(xasix~yasix),col=2,lty=1)

# exp3
print("START EXP3")
res <- associations_test(O, y, output$refactor_components[,1:K])
plot(res$x, res$y, main="refactor", xlab="uniform distribution", ylab="observed pvalues", pch='.', xlim=c(0,ceiling(max(res$x))),  ylim=c(0,ceiling(max(res$y))))
abline(lm(xasix~yasix),col=2,lty=1)

# exp4
print("START EXP4")
res <- associations_test(O, y, output$first_pca);
plot(res$x, res$y, main="P", xlab="uniform distribution", ylab="observed pvalues", pch='.', xlim=c(0,ceiling(max(res$x))),  ylim=c(0,ceiling(max(res$y))))
abline(lm(xasix~yasix),col=2,lty=1)

dev.off()
print("plot saved to plot.png file") #TODO move plot.png to argument
                                     #TODO change plot.png by time and date
system("open plot.png")

