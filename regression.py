# Y is an n*1 vector (numpy's ndarray)
# X is an n*d matrix (numpy's ndarray)

from sklearn import linear_model
from numpy import column_stack, zeros    

class LinearRegression(object):
    def __init__(self, y, x):
        zeros_vector =  zeros((len(x),  1))
    
        newx = column_stack((zeros_vector, x))

        regr = linear_model.LinearRegression()
        model = regr.fit(newx, y)
        pred = regr.predict(newx)
        self.residuals = y - pred
        self.coef = regr.coef_

class LogisticRegression(object):
    pass