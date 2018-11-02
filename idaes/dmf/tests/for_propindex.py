from idaes.core.property_base import PropertyParameterBase


class IndexMePlease1(PropertyParameterBase):

    @classmethod
    def define_metadata(cls, m):
        m.add_default_units({'temperature': 'K'})
        m.add_properties({'pressure': {'units': 'Pa', 'method': 'foo'},
                          'temperature': {'method': 'bar'}})
