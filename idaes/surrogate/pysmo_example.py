# from idaes.surrogate 
import pysmo_surrogate as ps
import pandas as pd
from io import StringIO
from idaes.surrogate.surrogate_block import SurrogateBlock


# df_test = pd.read_csv(r'C:\Users\OOAmusat\Downloads\rubbish_data.csv', header=0)

training_data = {'x1': [1, 2, 3, 4, 5], 'x2': [5, 6, 7, 8, 9], 'z1': [10, 20, 30, 40, 50], 'z2': [6, 8, 10, 12, 14]}
training_data = pd.DataFrame(training_data)
validation_data = {'x1': [1, 2, 3, 4], 'x2': [5, 6, 7, 8], 'z1': [10, 20, 30, 40], 'z2': [6, 8, 10, 12]}#{'x1': [2.5], 'x2': [6.5], 'z1': [25], 'z2': [9]}
validation_data = pd.DataFrame(validation_data)
input_labels = ["x1", "x2"]
output_labels = ["z1", "z2"]
bnds = {"x1": (0, 5), "x2": (0, 10)}

pysmo_trainer = ps.PysmoPolyTrainer(
    input_labels=input_labels,
    output_labels=output_labels,
    input_bounds=bnds,
    training_dataframe=training_data,
    validation_dataframe=validation_data,
    maximum_polynomial_order = 1,
    multinomials=False,
    number_of_crossvalidations=3,
    # extra_features = ['log(x1)', 'sin(x2)'],
    extra_features = ['x1 / x2'],
     # solution_method='pyomo'
    )

# Train surrogate
a1 = pysmo_trainer.train_surrogate()

p1 = ps.PysmoSurrogate(a1, input_labels, output_labels)
# Test surrogate evaluation
b1 = p1.evaluate_surrogate(validation_data[['x1', 'x2']])
# Test block generation
blk1 = SurrogateBlock(concrete=True)
blk1.build_model(p1)
blk1.pprint()
# Test save and load
strm = StringIO()
zv1 = p1.save(strm)
# strm.seek(0)
p1_load = ps.PysmoSurrogate.load(strm) # p1_load = p1.load(zv1)
c1 = p1_load.evaluate_surrogate(validation_data[['x1', 'x2']])
assert b1.equals(c1)
blk1_load = SurrogateBlock(concrete=True)
blk1_load.build_model(p1_load)
blk1_load.pprint()



pysmo_trainer = ps.PysmoRBFTrainer(
    input_labels=input_labels,
    output_labels=output_labels,
    input_bounds=bnds,
    training_dataframe=training_data,
    validation_dataframe=validation_data,
    solution_method = 'bfgs',
    regularization=False,
    basis_function = 'cubic')

# Train surrogate
a2 = pysmo_trainer.train_surrogate()

p2 = ps.PysmoSurrogate(a2, input_labels, output_labels, bnds)
# Test surrogate evaluation
b2 = p2.evaluate_surrogate(validation_data[['x1', 'x2']])
# Test block generation
blk2 = SurrogateBlock(concrete=True)
blk2.build_model(p2)
blk2.pprint()
# Test save and load
strm = StringIO()
zv2 = p2.save(strm)
# strm.seek(0)
p2_load = ps.PysmoSurrogate.load(strm)# p2_load = p2.load(zv2)
c2 = p2_load.evaluate_surrogate(validation_data[['x1', 'x2']])
assert b2.equals(c2)
blk2_load = SurrogateBlock(concrete=True)
blk2_load.build_model(p2_load)
blk2_load.pprint()



pysmo_trainer = ps.PysmoKrigingTrainer(
    input_labels=input_labels,
    output_labels=output_labels,
    numerical_gradients = True,
    regularization=False,
    # input_bounds=bnds,
    training_dataframe=training_data,
    validation_dataframe=validation_data)

# Train surrogate
a3 = pysmo_trainer.train_surrogate()

p3 = ps.PysmoSurrogate(a3, input_labels, output_labels, bnds)
# Test surrogate evaluation
b3 = p3.evaluate_surrogate(validation_data[['x1', 'x2']])
# Test block generation
blk3 = SurrogateBlock(concrete=True)
blk3.build_model(p3)
blk3.pprint()
# Test save and load
strm = StringIO()
zv3 = p3.save(strm)
# strm.seek(0)
p3_load = ps.PysmoSurrogate.load(strm) # p3_load = p3.load(zv3)
c3 = p3_load.evaluate_surrogate(validation_data[['x1', 'x2']])
assert b3.equals(c3)
blk3_load = SurrogateBlock(concrete=True)
blk3_load.build_model(p3_load)
blk3_load.pprint()

