from sklearn import preprocessing
from sklearn.decomposition import PCA as pca

class PCA( object ):
    def __init__(self, A):
        scaled = preprocessing.StandardScaler().fit(A).transform(A)
        pca_res = pca().fit(scaled) 

        self.U = pca_res.components_.transpose() # loadings
        self.P = pca_res.transform( scaled )     # scores

