#################################################
m = pyo.ConcreteModel()
m.TIME = pyo.Set(initialize=np.arange(1,5), ordered=True)
data_args = {1:[30.0],2:[20.0],3:[15.0],4:[30.0]}
def block_rule(b,t):
    b.process = create_mp_rankine_block(*data_args[t])
m.blocks = pyo.Block(m.TIME, rule=block_rule)


m.TIME.add(5)
m.blocks.keys().append(5)
m.blocks.values().append(create_mp_rankine_block(30.0))
m.blocks[5].transfer_attributes_from(create_mp_rankine_block(30.0))

def get_active_process_blocks(self):
    return [b.process for b in self.pyomo_model.blocks.values() if b.process.active]
    #return [p for p in self.pyomo_model.blocks[:].process if p.active]
    #return [self.pyomo_model.blocks[i].process for i in self.pyomo_model.TIME if self.pyomo_model.blocks[i].process.active]

m = pyo.ConcreteModel()
m.x = pyo.Var()
m.y = pyo.Var()
m.test = pyo.Constraint(range(2))

    # #fix variables on first block in the horizon (i.e. the current time)
    # def _fix_initial_states(self,b1,variable_pairs):
    #     b1.fix_states = pyo.Constraint(range(len(variable_pairs)))
    #     for (i,pair) in enumerate(variable_pairs):
    #         b1.fix_states[i] = pair[1]==pyo.value(pair[0])
