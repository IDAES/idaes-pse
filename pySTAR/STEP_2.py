#Generate a group of candidate expression trees using Yr as probability
#Starting from the root node
#As a parent/sibling node is determined the probability of the rest of the nodes change

from math import log2
from pyomo.environ import *
from copy import deepcopy
import numpy as np
import re

np.random.seed(42) 

def candidate_trees(X, mr, Yr, Br, Ur, Lr, NnotTr, Tr, Nr, num_trees):

    #Keep hyperparameters as same as STEP_1
    depth = log2(max(Tr)+1)-1
    Or = Br + Ur + Lr    
    Nodes = Nr

    #Converting Yr into probabilities
    Options_for_nodes, Ys_for_nodes = OptionsForNodes(mr, Nodes, Yr, Or)


    ### Tree Generation ##########################################################
    regexp = re.compile("cst[\d]*")
        
    rr_results = []
    for tree in range(num_trees):

        probs_for_nodes = deepcopy(Ys_for_nodes)

        VARS_with_0_values = [
            "x_" + str(i + 1) for i in range(X.shape[1]) if (X.T[i] == 0).any() #if (X.T[i] < 1e-2).any() and (X.T[i] >= 0).any()
        ]
        VARS_with_Negative_values = [
            "x_" + str(i + 1) for i in range(X.shape[1]) if (X.T[i] < 0).any()
        ]

        Yfix = {}
        def deleteChildren(n):
            if 2 * n in Nodes:
                Yfix[2 * n] = 0
                Yfix[2 * n + 1] = 0
                deleteChildren(2 * n)
                deleteChildren(2 * n + 1)
        
        # Yfix = {}
        # def descendents_nonzero(n):
        #     for  in Nodes:
                

        for node in Nodes:
            
            l = 2 * node
            r = 2 * node + 1
            if node > 1:
                if node % 2 == 0:
                    parent = int(node / 2)
                else:
                    parent = int((node - 1) / 2)
            else:
                parent = None

            if node not in Yfix.keys():
                options = Options_for_nodes[node]

                try:
                    choice = np.random.choice(
                        np.arange(len(options)), p=probs_for_nodes[node]
                    )
                except:
                    pass

                node_value_chosen = options[choice]
                Yfix[node] = node_value_chosen
                value_chosen = node_value_chosen[1]

                if value_chosen == "/":
                    if len(VARS_with_0_values) == X.shape[1] and r in Tr:
                        for idx, o in enumerate(Options_for_nodes[node]):
                            if o[1] == "/":
                                
                                probs_for_nodes[node][idx] = 0
                                if np.sum(probs_for_nodes[node]) == 0:
                                    probs_for_nodes[node] = np.zeros(len(probs_for_nodes[node]))
                                else: probs_for_nodes[node] = probs_for_nodes[node] / np.sum(probs_for_nodes[node])
                                
                                choice = np.random.choice(
                                    np.arange(len(options)), p=probs_for_nodes[node]
                                )
                                node_value_chosen = options[choice]
                                Yfix[node] = node_value_chosen
                                value_chosen = node_value_chosen[1]
                    else:
                        for idx, o in enumerate(Options_for_nodes[r]):
                            if o[1] == "cst" or o[1] in VARS_with_0_values:
                                probs_for_nodes[r][idx] = 0
                                if np.sum(probs_for_nodes[r]) == 0:
                                    probs_for_nodes[r] = np.zeros(len(probs_for_nodes[r]))
                                else: probs_for_nodes[r] = probs_for_nodes[r] / np.sum(probs_for_nodes[r])
                    # else:
                    #     for idx, o in enumerate(Options_for_nodes[r]):
                    #         if o[1] in VARS_with_0_values:
                                
                                
                if value_chosen == "**0.5":
                    if (
                        #len(np.union1d(VARS_with_Negative_values, VARS_with_0_values))
                        len(VARS_with_Negative_values) == X.shape[1]
                    ):
                        for idx, o in enumerate(options):
                            if o[1] == "**0.5":
                                probs_for_nodes[node][idx] = 0
                                if np.sum(probs_for_nodes[node]) == 0:
                                    probs_for_nodes[node] = np.zeros(len(probs_for_nodes[node]))
                                else: probs_for_nodes[node] = probs_for_nodes[node] / np.sum(probs_for_nodes[node])
                                
                                choice = np.random.choice(
                                    np.arange(len(options)), p=probs_for_nodes[node]
                                )
                                node_value_chosen = options[choice]
                                Yfix[node] = node_value_chosen
                                value_chosen = node_value_chosen[1]
                    else:
                        for idx, o in enumerate(Options_for_nodes[r]):
                            if (
                                o[1] == "cst" or o[1] in VARS_with_Negative_values
                                #or o[1] in VARS_with_Negative_values + VARS_with_0_values
                            ):
                                probs_for_nodes[r][idx] = 0
                                if np.sum(probs_for_nodes[r]) == 0:
                                    probs_for_nodes[r] = np.zeros(len(probs_for_nodes[r]))
                                else: probs_for_nodes[r] = probs_for_nodes[r] / np.sum(probs_for_nodes[r])
                                
                if value_chosen == "log":
                    if (
                        len(np.union1d(VARS_with_Negative_values, VARS_with_0_values))
                        == X.shape[1]
                    ):
                        for idx, o in enumerate(options):
                            if o[1] == "log":
                                probs_for_nodes[node][idx] = 0
                                if np.sum(probs_for_nodes[node]) == 0:
                                    probs_for_nodes[node] = np.zeros(len(probs_for_nodes[node]))
                                else: probs_for_nodes[node] = probs_for_nodes[node] / np.sum(probs_for_nodes[node])
                                
                                choice = np.random.choice(
                                    np.arange(len(options)), p=probs_for_nodes[node]
                                )
                                node_value_chosen = options[choice]
                                Yfix[node] = node_value_chosen
                                value_chosen = node_value_chosen[1]
                    else:
                        for idx, o in enumerate(Options_for_nodes[r]):
                            if (
                                o[1] == "cst"
                                or o[1] in VARS_with_Negative_values + VARS_with_0_values
                            ):
                                probs_for_nodes[r][idx] = 0
                                if np.sum(probs_for_nodes[r]) == 0:
                                    probs_for_nodes[r] = np.zeros(len(probs_for_nodes[r]))
                                else: probs_for_nodes[r] = probs_for_nodes[r] / np.sum(probs_for_nodes[r])
                                
                if value_chosen not in ["*", "+"] and node not in Tr:
                    for idx, o in enumerate(Options_for_nodes[r]):
                        if o[1] == "cst":
                            probs_for_nodes[r][idx] = 0
                            if np.sum(probs_for_nodes[r]) == 0:
                                probs_for_nodes[r] = np.zeros(len(probs_for_nodes[r]))
                            else: probs_for_nodes[r] = probs_for_nodes[r] / np.sum(probs_for_nodes[r])
                
                if value_chosen == "-":
                    for idx, o in enumerate(Options_for_nodes[r]):
                        if o[1] == "cst":
                            probs_for_nodes[r][idx] = 0
                            if np.sum(probs_for_nodes[r]) == 0:
                                probs_for_nodes[r] = np.zeros(len(probs_for_nodes[r]))
                            else: probs_for_nodes[r] = probs_for_nodes[r] / np.sum(probs_for_nodes[r])
                            
                    for idx, o in enumerate(Options_for_nodes[r]):
                        if o[1] == "cst":
                            probs_for_nodes[r][idx] = 0
                            if np.sum(probs_for_nodes[r]) == 0:
                                probs_for_nodes[r] = np.zeros(len(probs_for_nodes[r]))
                            else: probs_for_nodes[r] = probs_for_nodes[r] / np.sum(probs_for_nodes[r])
                            
                if value_chosen in Ur:
                    for idx, o in enumerate(Options_for_nodes[r]):
                        if o[1] == "cst":
                            probs_for_nodes[r][idx] = 0
                            if np.sum(probs_for_nodes[r]) == 0:
                                probs_for_nodes[r] = np.zeros(len(probs_for_nodes[r]))
                            else: probs_for_nodes[r] = probs_for_nodes[r] / np.sum(probs_for_nodes[r])
                    Yfix[l] = 0
                    deleteChildren(l)
                
                if value_chosen == "cst" and node % 2 == 0:
                    for idx, o in enumerate(Options_for_nodes[node + 1]):
                        if o[1] == "cst":
                            probs_for_nodes[node + 1][idx] = 0
                            if np.sum(probs_for_nodes[node+1]) == 0:
                                probs_for_nodes[node+1] = np.zeros(len(probs_for_nodes[node+1]))
                            else: probs_for_nodes[node+1] = probs_for_nodes[node+1] / np.sum(probs_for_nodes[node+1])
                            
                if value_chosen in Lr and node not in Tr:
                    deleteChildren(node)
                
                if value_chosen in Lr and node % 2 == 0 and Yfix[parent][1] == "-":
                    for idx, o in enumerate(Options_for_nodes[node + 1]):
                        if o[1] == value_chosen:
                            probs_for_nodes[node + 1][idx] = 0
                            if np.sum(probs_for_nodes[node+1]) == 0:
                                probs_for_nodes[node+1] = np.zeros(len(probs_for_nodes[node+1]))
                            else: probs_for_nodes[node+1] = probs_for_nodes[node+1] / np.sum(probs_for_nodes[node+1])
                            

        es, cst_count = YfixToExpression(regexp, Yfix, depth, enhance_subtree=True, enhance_var=True)

        rr_results.append((Yfix, es, cst_count))

    return rr_results


def OptionsForNodes(mr, Nodes, Yr, Or):
    
    #initialize
    Options_for_nodes, Ys_for_nodes = (
        {n: {o: False for o in Or} for n in Nodes},
        {n: {o: 0 for o in Or} for n in Nodes},
    )

    for node in Nodes:
        Options_for_node = [k for k in Yr if k[0] == node]
        Options_for_nodes[node] = Options_for_node
        
        if np.sum([mr.y[k].value for k in Options_for_node]) == 0:
            Ys_for_node = [0]
        else:
            Ys_for_node = [mr.y[k].value / np.sum([mr.y[k].value for k in Options_for_node]) for k in Options_for_node]
        Ys_for_nodes[node] = Ys_for_node
        
    return Options_for_nodes, Ys_for_nodes

def YfixToExpression(regexp, Yfix, depth, enhance_var=True, enhance_subtree=True):
    
    Yfix_tree = {}
    for k in Yfix.keys():
        if Yfix[k] != 0:
            Yfix_tree[k] = Yfix[k][1]

    es = TreeToExpressionString(
        Yfix_tree,
        node=1,
        depth=depth,
        enhance_var=enhance_var,
        enhance_subtree=enhance_subtree,
    )
    cst_count = len(regexp.findall(es))
    return es, cst_count

def TreeToExpressionString(tree, node, depth, enhance_var=True, enhance_subtree=True):
    if node * 2 in tree.keys() and node * 2 + 1 in tree.keys():
        
        l = TreeToExpressionString(tree, node * 2, depth, enhance_var, enhance_subtree)
        r = TreeToExpressionString(tree, node * 2 + 1, depth, enhance_var, enhance_subtree)
        
        p = tree[node]

        current = "+"
        
        if p == "+":
            current = "+"
        elif p == "-":
            current = "-"
        elif p == "*":
            current = "*"
        elif p == "/":
            current = "/"
        elif p == "**0.5":
            current = "**0.5"
        #         elif p == "sqr":
        #             current = '**2'
        #         elif p == "cub":
        #             current = '**3'
        if enhance_subtree:
            return "(cst + cst*(%s%s%s))" % (l, current, r)
        else:
            return "(%s%s%s)" % (l, current, r)
    elif node * 2 + 1 in tree.keys():
        l = ""
        r = TreeToExpressionString(
            tree, node * 2 + 1, depth, enhance_var, enhance_subtree
        )
        p = tree[node]

        current = "+"
        if p == "**0.5":
            current = "**0.5"
            if enhance_subtree:
                return "(cst + cst*(%s%s)%s)" % (l, r, current)
            else:
                return "(%s%s)%s" % (l, r, current)
        #         elif p == "sqr":
        #             current = '**2'
        #             return "(%s%s)%s" % (l, r, current)
        #         elif p == "cub":
        #             current = '**3'
        #             return "(%s%s)%s" % (l, r, current)
        elif p == "log":
            current = "log("
            if enhance_subtree:
                return "cst*%s%s%s)" % (current, l, r)
            else:
                return "(%s%s%s))" % (current, l, r)
        elif p == "exp":
            current = "exp("
            if enhance_subtree:
                return "cst*%s%s%s)" % (current, l, r)
            else:
                return "%s%s%s)" % (current, l, r)
        if enhance_subtree:
            return "(cst*(%s%s%s))" % (l, current, r)
        else:
            return "(%s%s%s)" % (l, current, r)
    else:
        if node in tree.keys():
            if tree[node][0] == "x":
                if enhance_var or enhance_subtree:
                    if node % 2 == 0:  # left node
                        return "cst + cst*X[:," + tree[node][2:] + "-1]"
                    else:
                        return "(cst*X[:," + tree[node][2:] + "-1]) + cst"
                else:
                    return "X[:," + tree[node][2:] + "-1]"
            else:
                return tree[node]


