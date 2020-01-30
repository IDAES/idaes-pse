# MatOpt

MatOpt is a package for facilitating the design of nanostructured materials via mathematical optimization.
We provide features that can be grouped into two main contributions. 
First, we provide several objects for representing the nanostructured materials design space and for constructing data structures useful in generic optimization. 
Second, we provide a modeling framework for specifying material optimization problems and casting them as mixed integer linear programming optimization models. 

For details of the methodology, see PAPER. 

## Examples

In the examples folder, we provide five case studies that can be modeled and solved by MatOpt.
In each case, we provide a Jupyter notebook with explanation as well as an equivalent Python script.

## Dependencies

1. Pyomo
   We use Pyomo for constructing the optimization models. 
2. Numpy 
   We use Numpy arrays to represent geometric points.
3. CPLEX
   The tool has been optimized to work with the MIP solver CPLEX. Users who do not have access to CPLEX via Pyomo are advised to use NEOS-CPLEX. This can be selected by edditing the source code for MatOptModel.
