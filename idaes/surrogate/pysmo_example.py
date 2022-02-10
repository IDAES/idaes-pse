# from idaes.surrogate 
import pysmo_surrogate as ps
import pandas as pd
from io import StringIO
from idaes.surrogate.surrogate_block import SurrogateBlock


# df_test = pd.read_csv(r'C:\Users\OOAmusat\Downloads\rubbish_data.csv', header=0)

training_data = {'x1': [1, 2, 3, 4], 'x2': [5, 6, 7, 8], 'z1': [10, 20, 30, 40], 'z2': [6, 8, 10, 12]}
training_data = pd.DataFrame(training_data)
validation_data = {'x1': [1, 2, 3, 4], 'x2': [5, 6, 7, 8], 'z1': [10, 20, 30, 40], 'z2': [6, 8, 10, 12]}#{'x1': [2.5], 'x2': [6.5], 'z1': [25], 'z2': [9]}
validation_data = pd.DataFrame(validation_data)
input_labels = ["x1", "x2"]
output_labels = ["z1", "z2"]
bnds = {"x1": (0, 5), "x2": (0, 10)}

pysmo_trainer = ps.PysmoPolyTrainer(
    input_labels=input_labels,
    output_labels=output_labels,
    # input_bounds=bnds,
    training_dataframe=training_data,
    validation_dataframe=validation_data,
    maximum_polynomial_order = 1, 
    multinomials=False,
    number_of_crossvalidations=5,
    solution_method='pyomo')

a = pysmo_trainer.train_surrogate()


p2 = ps.PysmoSurrogate(a, input_labels, output_labels)
# Test surrogate evaluation
b = p2.evaluate_surrogate(validation_data[['x1', 'x2']])
# Test block generation
blk = SurrogateBlock(concrete=True)
blk.build_model(p2)

strm = StringIO()
p2.save(strm)
print(strm)
ps.PysmoSurrogate.load(strm)



pysmo_trainer = ps.PysmoRBFTrainer(
    input_labels=input_labels,
    output_labels=output_labels,
    # input_bounds=bnds,
    training_dataframe=training_data,
    validation_dataframe=validation_data,
    solution_method = 'bfgs',
    regularization=False,
    basis_function = 'cubic')

a = pysmo_trainer.train_surrogate()


p2 = ps.PysmoSurrogate(a, input_labels, output_labels)
# Test surrogate evaluation
b = p2.evaluate_surrogate(validation_data[['x1', 'x2']])
# Test block generation
blk = SurrogateBlock(concrete=True)
blk.build_model(p2)




pysmo_trainer = ps.PysmoKrigingTrainer(
    input_labels=input_labels,
    output_labels=output_labels,
    # input_bounds=bnds,
    training_dataframe=training_data,
    validation_dataframe=validation_data)

a = pysmo_trainer.train_surrogate()


p2 = ps.PysmoSurrogate(a, input_labels, output_labels)
# Test surrogate evaluation
b = p2.evaluate_surrogate(validation_data[['x1', 'x2']])
# Test block generation
blk = SurrogateBlock(concrete=True)
blk.build_model(p2)

# import pickle
# import json
# import sys
# import os

# # open pickle file
# obj = open('z1_solution.pickle', 'rb')
            
# # convert pickle object to json object
# json_obj = json.loads(json.dumps(obj, default=str))
# print(json_obj)

# # write the json file
# with open(
#         os.path.splitext(sys.argv[1])[0] + '.json',
#         'w',
#         encoding='utf-8'
#     ) as outfile:
#     json.dump(json_obj, outfile, ensure_ascii=False, indent=4)
