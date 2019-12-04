
from abc import abstractmethod

class BBlock(object):   
    """ """
    # === PROPERTY EVALUATION METHODS
    @abstractmethod
    def __eq__(self,other): 
        raise NotImplementedError

    @abstractmethod
    def __le__(self,other):
        raise NotImplementedError



