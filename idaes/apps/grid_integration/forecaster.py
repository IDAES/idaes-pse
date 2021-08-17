import pandas as pd

class WhiteNoiseForecaster:

    def __init__(self, price_forecasts_df):
        self.price_forecasts_df = price_forecasts_df
        self.price_forecasts_df.set_index(['Date', 'Hour'], inplace = True)
        self._rename_columns()

    def _rename_columns(self):
        col_mapping = {"Scenario {}".format(i): i - 1 for i in range(1, 11)}
        self.price_forecasts_df.rename(columns = col_mapping, inplace = True)

        return

    def forecast(self, date, **kwargs):
        date = str(date)
        return self.price_forecasts_df.loc[date].to_dict('list')


if __name__ == "__main__":

    price_forecasts_df = pd.read_csv('examples/lmp_forecasts_example.csv')
    forecaster = WhiteNoiseForecaster(price_forecasts_df = price_forecasts_df)

    date = "2021-08-01"
    forecasts = forecaster.forecast(date = date)
