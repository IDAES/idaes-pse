class FlowsheetNotFound(Exception):
    def __init__(self, id_, location):
        super().__init__(f"Flowsheet {id_} not found in {location}")
        self.location = location  # to help distinguish


class FlowsheetNotFoundInDatastore(FlowsheetNotFound):
    def __init__(self, id_):
        super().__init__(id_, "datastore")


class FlowsheetNotFoundInMemory(FlowsheetNotFound):
    def __init__(self, id_):
        super().__init__(id_, "Python process memory")


class FlowsheetUnknown(Exception):
    def __init__(self, id_):
        super().__init__(f"Unrecognized flowsheet '{id_}'")


class ProcessingError(Exception):
    """Use for errors processing input."""
    pass


class VisualizerError(Exception):
    """Generic error for visualizer
    """
    pass


class VisualizerSaveError(VisualizerError):
    def __init__(self, save_as, error_message):
        msg = f"While saving flowsheet as '{save_as}': {error_message}"
        super().__init__(msg)


class DatastoreError(Exception):
    pass


class DatastoreSaveError(DatastoreError):
    pass