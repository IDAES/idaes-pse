class TrainingMetrics(object):
    def __init__(self, RMSE, MSE, SSE, R2):
        # root mean-squared error
        self.RMSE = RMSE

        # mean squared error
        self.MSE = MSE

        # sum of squared errors
        self.SSE = SSE

        # R-squared
        self.R2 = R2

        # TODO: number of datapoints out of bounds
 
        
# TODO: Functions or methods to get metrics from training and validation data
#def compute_fit_metrics(surrogate_object, data):
#    # do the stuff
#    return TrainingMetrics(stuff)

