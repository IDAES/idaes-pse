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


class DatastoreSerializeError(DatastoreError):
    def __init__(self, obj, err, stream=None):
        to_stream = "" if stream is None else f" to '{stream}'"
        message = f"Serializing object {obj}{to_stream} as JSON failed: {err}"
        super().__init__(message)


class DatastoreSaveError(DatastoreError):
    pass


class TooManySavedVersions(Exception):
    pass