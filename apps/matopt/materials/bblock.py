from abc import abstractmethod


class BBlock(object):
    """An abstract class for material building blocks."""

    # === PROPERTY EVALUATION METHODS
    @abstractmethod
    def __eq__(self, other):
        raise NotImplementedError

    @abstractmethod
    def __le__(self, other):
        raise NotImplementedError
