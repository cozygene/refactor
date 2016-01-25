from numpy import loadtxt, delete, isnan, nanvar, where, column_stack, ones
from sklearn import linear_model

class Regress( object ):
	
	def __init__(self):
		return None
	
	@staticmethod
	def regress(y,x):
		"""
		Regresses y on x and returns the residuals.
		"""
		ones_vector =  ones(len(y))
		newx = column_stack((ones_vector, x))
		regr = linear_model.LinearRegression()
		model = regr.fit(newx, y)
		pred = regr.predict(newx)
		res = (y - pred)
		return res

