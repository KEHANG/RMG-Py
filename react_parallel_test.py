from rmgpy.rmg.main import RMG
from rmgpy.species import Species
from rmg import parseCommandLineArguments

if __name__ == '__main__':

    # parse the args
    args = parseCommandLineArguments()
    # create RMG instance
    rmg = RMG()
    # initialize the rmg job
    # which includes loading inputfile
    # loading database and core&edge model
    rmg.initialize(args)

    # load new core species
    newCoreSpecies = Species().fromSMILES('CC(C)C')

    # enlarge model by running react
    # running react in serial mode or parallel mode
    rmg.reactionModel.enlarge(newCoreSpecies)

    # save everything
    rmg.saveEverything()
