from rmgpy.rmg.main import *

def loadRestartFile(rmg, path):
    """
    Load a restart file at `path` on disk.
    """

    import cPickle

    # Unpickle the reaction model from the specified restart file
    print('Loading previous restart file...')
    f = open(path, 'rb')
    rmg.reactionModel = cPickle.load(f)
    f.close()

    # A few things still point to the species in the input file, so update
    # those to point to the equivalent species loaded from the restart file
    # The reactions and reactionDict still point to the old reaction families
    reactionDict = {}
    for family0 in rmg.reactionModel.reactionDict:

        # Find the equivalent library or family in the newly-loaded kinetics database
        family = None
        if isinstance(family0, KineticsLibrary):
            for label, database in rmg.database.kinetics.libraries.iteritems():
                if database.label == family0.label:
                    family = database
                    break
        elif isinstance(family0, KineticsFamily):
            for label, database in rmg.database.kinetics.families.iteritems():
                if database.label == family0.label:
                    family = database
                    break
        else:
            import pdb; pdb.set_trace()
        if family is None:
            raise Exception("Unable to find matching reaction family for %s" % family0.label)

        # Update each affected reaction to point to that new family
        # Also use that new family in a duplicate reactionDict
        reactionDict[family] = {}
        for reactant1 in rmg.reactionModel.reactionDict[family0]:
            reactionDict[family][reactant1] = {}
            for reactant2 in rmg.reactionModel.reactionDict[family0][reactant1]:
                reactionDict[family][reactant1][reactant2] = []
                if isinstance(family0, KineticsLibrary):
                    for rxn in rmg.reactionModel.reactionDict[family0][reactant1][reactant2]:
                        assert isinstance(rxn, LibraryReaction)
                        rxn.library = family
                        reactionDict[family][reactant1][reactant2].append(rxn)
                elif isinstance(family0, KineticsFamily):
                    for rxn in rmg.reactionModel.reactionDict[family0][reactant1][reactant2]:
                        assert isinstance(rxn, TemplateReaction)
                        rxn.family = family
                        reactionDict[family][reactant1][reactant2].append(rxn)

    rmg.reactionModel.reactionDict = reactionDict

# the purpose of this script is to line-profile processNewReactions
if __name__ == '__main__':
    # to run processNewReactions we need 3 inputs
    # newReactions, newSpecies, pdepnetworks (which can be None here)

    # set-up RMG object
    rmg = RMG()
    # load kinetic database and forbidden structures
    rmg.database = RMGDatabase()
    path = os.path.join(os.path.dirname(__file__), '..', 'RMG-database', 'input')
    print 'Path is:', path
    print("loading forbidden structures...")
    rmg.database.loadForbiddenStructures(os.path.join(path, 'forbiddenStructures.py'))
    print("loading kinetics families...")
    rmg.database.loadKinetics(os.path.join(path, 'kinetics'), kineticsFamilies='default')
    print "Succeeded in loading database!"

    # unpickle restart.pkl
    restartFilePath = os.path.join(os.getcwd(), 'restart110.pkl')
    loadRestartFile(rmg, restartFilePath)
    
    # newSpecies
    for species in rmg.reactionModel.edge.species:
        if species.index == 4691:
            newSpecies = species
            break
    
    # react
    print "Start react!"
    newReactions = []
    newReactions.extend(rmg.reactionModel.react(rmg.database, newSpecies))
    for coreSpecies in rmg.reactionModel.core.species:
        if coreSpecies.reactive:
            newReactions.extend(rmg.reactionModel.react(rmg.database, newSpecies, coreSpecies))
    newReactions.extend(rmg.reactionModel.react(rmg.database, newSpecies, newSpecies))
    print "react finished!"
    # processNewReaction
    # Add new species
    reactionsMovedFromEdge = rmg.reactionModel.addSpeciesToCore(newSpecies)
    
    # Process the new reactions
    # While adding to core/edge/pdep network, this clears atom labels:
    rmg.reactionModel.processNewReactions(newReactions, newSpecies, None)

    corerxns = rmg.reactionModel.core.reactions
    edgerxns = rmg.reactionModel.edge.reactions
    edgeSpe = rmg.reactionModel.edge.species
    print len(corerxns), len(edgerxns), len(edgeSpe)
    
