
if(!exists("refactor", mode="function")) source("refactor.R")

associations_test <- function(O, y, model_append)
{
    observed_pvalues <- c()
    for (site in 1:ncol(O))
    {
        model <- lm(y ~ O[,site] + model_append)
        
        pvalue <- coef(summary(model))[2,4]
        observed_pvalues[site] = as.numeric(pvalue)
    }

    return(observed_pvalues) 
}

draw_qqplot <- function(y, title, xtitle, ytitle, style='.')
{
    x <- runif(length(y))
    x <- sort(-log10(x))
    y <- sort(-log10(y)) 
    plot(x, y, main=title, xlab=xtitle, ylab=ytitle, pch=style, xlim=c(0,ceiling(max(x))),  ylim=c(0,ceiling(max(y))))
    
    # add y=x trend
    xasix<-0:10
    yasix<-0:10
    abline(lm(xasix~yasix),col=2,lty=1)
}

args <- commandArgs(trailingOnly = TRUE)
data_file=args[1] # data file
y=args[2] # a phenotype file
K=as.numeric(args[3]) # an integer
R=args[4] # a matrix with k columns; these are the cell proportions)


#TODO enter to funs?
    O = as.matrix(read.table(data_file))    
    sample_id_O <- O[1, -1] # extract samples ID
    O <- O[-1,] # remove sample ID from matrix
    cpgnames <- O[, 1] ## set rownames
    O <- O[, -1] 
    O = t(matrix(as.numeric(O),nrow=nrow(O),ncol=ncol(O)))

phenotype_matrix = as.matrix(read.table(y))
y <-  matrix(as.numeric(as.matrix(phenotype_matrix[, -1]) ))



png('plot.png')
par(mfrow=c(2,2))


output <- refactor(data_file, K)
print("Unadjusted analysis...")
# exp1
print("START EXP1")
observed_pvalues <- associations_test(O,y, matrix(0, nrow = nrow(O), ncol = 1))
draw_qqplot(observed_pvalues, title='Unadjusted analysis', xtitle='-log10(expected)', ytitle='-log10(observed)')

# exp2
print("START EXP2")
observed_pvalues <- associations_test(O, y, as.matrix(read.table(R)))
draw_qqplot(observed_pvalues, title='Adjusted analysis using cell proportions', xtitle='-log10(expected)', ytitle='-log10(observed)')

# exp3
print("START EXP3")
observed_pvalues <- associations_test(O, y, output$refactor_components[,1:K])
draw_qqplot(observed_pvalues, title='Adjusted analysis using ReFACTor', xtitle='-log10(expected)', ytitle='-log10(observed)')

# exp4
print("START EXP4")
observed_pvalues <- associations_test(O, y, output$first_pca);
draw_qqplot(observed_pvalues, title='Adjusted analysis using PCA', xtitle='-log10(expected)', ytitle='-log10(observed)')

dev.off()
print("plot saved to plot.png file") #TODO move plot.png to argument
                                     #TODO change plot.png by time and date
system("open plot.png")

