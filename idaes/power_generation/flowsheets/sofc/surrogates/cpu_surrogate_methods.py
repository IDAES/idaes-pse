##############################################################################
# The development of this flowsheet/code is funded by the ARPA-E DIFFERENTIATE
# project: “Machine Learning for Natural Gas to Electric Power System Design”
# Project number: DE-FOA-0002107-1625.
# This project is a collaborative effort between the Pacific Northwest National
# Laboratory, the National Energy Technology Laboratory, and the University of
# Washington to design NGFC systems with high efficiencies and low CO2
# emissions.
##############################################################################


class compressor_power_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 16592.309504827182536246 * x1 + 17806584.568199951201677 * x3 - 15529374.443028066307306


class heat_duty_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 27431.344430845900205895 * x1 + 25040167.356691692024469 * x3 - 21837134.557537198066711


class refrigeration_duty_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return - 0.30738967679398379205177E-004 * x5 + 0.45666609499185339089765E-006


class pureco2_flow_mol_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.84491492199656215156267 * x1 + 1381.8486099781443954271 * x3 - 1204.9518833417739642755


class pureco2_ar_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.87539452852376016910194E-007 * x1 + 0.31019041198388018487631E-001 * x2 - 0.16384578403083416373726E-002 * x5 - 0.11401484796994311183421E-002 * x6


class pureco2_co2_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.84491226339466729555738 * x1 + 1381.8776912884593457420 * x3 - 1204.9772340007211823831


class pureco2_o2_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.25557869821781261075201E-005 * x1 - 0.30160307289038564698691E-001 * x3 + 0.26298816057506283622169E-001


class pureco2_h2o_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.00


class pureco2_n2_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.69246160597454559832190E-006 * x3 - 0.56005883473002498644045E-006


class pureco2_temperature_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 10.267585306350582641244 * x2 + 32.023906991174357017371 * x3 + 30.302544593098598824099 * x4 + 278.86088998572887476257


class pureco2_pressure_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return - 0.28666413419197245189076E-011 * x1 + 15271893.400000000372529


class water_flow_mol_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.64898275330398047522351E-001 * x1 + 995.23584715669267097837 * x4 - 64.595577093201725915605


class water_ar_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.23736607897806708417605E-008 * x1 + 0.83689507888515507049582E-003 * x2 - 0.14749984695782550662112E-003 * x5


class water_co2_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.36504473074529521385892E-005 * x1 + 0.39821667908625070844697E-002 * x3 + 0.30283921974618616818065E-001 * x4 - 0.54375872162499300915828E-002


class water_o2_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.46405004009352215417855E-007 * x1 - 0.99119387109367312720110E-003 * x3 + 0.86473506719933338534462E-003


class water_h2o_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.64894576227660510925332E-001 * x1 + 995.21157061827375400753 * x4 - 64.594007608678808196601


class water_n2_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.71432912211586127618749E-009 * x1


class water_temperature_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.47345180660722969989695E-002 * x3 - 1.1291949419088465056404 * x4 - 0.13703395454789617582958E-001 * x5 + 311.08988060641769379799


class water_pressure_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 117210.91999999999825377


class vent_flow_mol_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.90250925128294370858306E-001 * x1 - 1373.8223230605847220431 * x3 + 1197.7852419317096064333


class vent_ar_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.27178896806122057258626E-002 * x1 + 953.43145974385458885081 * x2 - 168.06049873202732669597 * x5


class vent_co2_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.27083917235717330090905E-001 * x1 - 407.85524106058562665567 * x3 + 355.59132055263916072363


class vent_o2_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.44846813052226472406936E-001 * x1 - 956.47536271014735120843 * x3 + 834.42603654916342748038


class vent_h2o_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.00


class vent_n2_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.15024528735008041077648E-001 * x1


class vent_temperature_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return - 0.43741475554195015242121E-016 * x1 + 305.37222222222220580079


class vent_pressure_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return - 0.22395635483747847803966E-013 * x1 + 115831.96800000000803266
