import pandas as pd
import numpy as np  
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('font', size=18)
plt.rc('axes', titlesize=18)
import os


# result_data_dir = 'prescient_verification_results/pmin_175_nn_free_startup_2/'
result_data_dir = 'bidding_plugin_test_multiperiod_rankine'


bus_detail_df = pd.read_csv(os.path.join(result_data_dir,'bus_detail.csv'))

# thermal detail (power dispatched by each generator)
thermal_detail_df = pd.read_csv(os.path.join(result_data_dir,'thermal_detail.csv'))

# renewable details
renewable_detail_df = pd.read_csv(os.path.join(result_data_dir,'renewables_detail.csv'))

# line detail (this has the power flow on each line)
line_detail_df = pd.read_csv(os.path.join(result_data_dir,'line_detail.csv'))

# hourly summary
hourly_summary_df = pd.read_csv(os.path.join(result_data_dir,'hourly_summary.csv'))

# the list of unique thermal generators
generator_list = pd.unique(thermal_detail_df['Generator'])

# the list of unique renewable power plants
renewable_list = pd.unique(renewable_detail_df['Generator'])

bidding_df = pd.read_csv(os.path.join(result_data_dir,'bidding_detail.csv'))
tracking_df = pd.read_csv(os.path.join(result_data_dir,'tracking_detail.csv'))

#Look at an example coal generator (This is the generator we studied for the big Prescient simulation data set)
coal_generator = '102_STEAM_3'
gen_results = thermal_detail_df.loc[thermal_detail_df['Generator'] == coal_generator] #output results for this generator
dispatch = gen_results["Dispatch"]


bid_time = bidding_df["Hour"]
bid_profiles = bidding_df.iloc[:,3:-1]
alphas = np.linspace(0,1,len(bid_time))


bid_profiles_sorted = bid_profiles.reindex(sorted(bid_profiles.columns), axis=1)

fig, ax = plt.subplots(figsize = (8,8))
ax.set_xlabel("Power Output [MW]", fontsize=24)
ax.set_ylabel("Cost [$]", fontsize=24)
for t in range(len(bid_time)):
    costs = bid_profiles_sorted.iloc[t,0:9].to_numpy()
    powers = bid_profiles_sorted.iloc[t,9:18].to_numpy()
    powers = powers[~np.isnan(powers)]
    costs = costs[~np.isnan(costs)]
    plt.plot(powers,costs,linewidth = 2, alpha = alphas[t],color = "black")

plt.tight_layout()
fig.savefig("bidding_results.png")
