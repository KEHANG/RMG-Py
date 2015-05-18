from rmgpy.rmg.main import *
import time
import sys
import json

if __name__ == '__main__':

    t0 = time.time()

    # set-up RMG object
    rmg = RMG()
    rmg.reactionModel = CoreEdgeReactionModel()

    # load kinetic database and forbidden structures
    rmg.database = RMGDatabase()
    path = os.path.join(os.path.dirname(__file__), '..', 'RMG-database', 'input')
    print 'Path is:', path
    print("loading forbidden structures...")
    rmg.database.loadForbiddenStructures(os.path.join(path, 'forbiddenStructures.py'))
    print("loading kinetics families...")
    rmg.database.loadKinetics(os.path.join(path, 'kinetics'), kineticsFamilies='default')

    # create a list of old core species from thermo libraries
    thermoLibrary = 'CHO'
    rmg.database.loadThermo(os.path.join(path, 'thermo'), thermoLibraries=[thermoLibrary])

    ## extract the entries out from database
    entryList = rmg.database.thermo.libraries[thermoLibrary].entries.values()
    entryNum = len(entryList)
    speciesList = []
    rootSpeciesDict = {} # for IDize part
    for i in range(entryNum):
        molecule = entryList[i].item
        index = entryList[i].index
        label = entryList[i].label
        species = Species(index = index, label = label, molecule = [molecule])
        speciesList.append(species)
        rootSpeciesDict[species.index] = species

    ## choose the 2nd species to final species to be oldCoreSpeciesList
    oldCoreSpeciesList = speciesList[1:]
    ### to make the old core species list large enough
    ### so that the execution time is long enough to overcome
    ### parallel overhead
    copyNum = int(sys.argv[2])
    for _ in range(copyNum-1):
        oldCoreSpeciesList.extend(speciesList[1:])

    ## choose the 1st species to be the new core species
    newSpecies = speciesList[0]

    # set-up step done
    t1 = time.time()
    print "Set-up step takes {0} secs.".format(t1-t0)

    # implement react method for newSpecies and oldCoreSpeciesList
    parallelMode = (sys.argv[1]=='parallel')
    newReactions = []
    entryStr = ''

    ## choose serial mode
    if not parallelMode:
        for coreSpecies in oldCoreSpeciesList:
            newReactions.extend(rmg.reactionModel.react(rmg.database, newSpecies, coreSpecies))
    ## choose parallel mode
    else:
        from scoop.futures import map
        from scoop import shared
        import gc

        ## broadcast database and forbiddenstructures
        families = rmg.database.kinetics.families
        print("Sharing kinetics families....")
        t_setConst_0 = time.time()
        shared.setConst(database_kinetics_families = families) # families is a dictionary
        print("Sharing forbiddenStructures...")
        shared.setConst(database_forbiddenStructures = rmg.database.forbiddenStructures)
        t_setConst_1 = time.time()
        t_setConst = t_setConst_1 - t_setConst_0
        print("Done broadcasting.")


        ## spawning tasks
        caseNum = 3
        corespeciesList = oldCoreSpeciesList
        corespeciesNum = len(corespeciesList)
        workerNum = int(sys.argv[3])
        taskNum = workerNum
        corespeciesList_list = []
        for i in range(taskNum):
            corespeciesList_list.append(corespeciesList[corespeciesNum*i/taskNum:corespeciesNum*(i+1)/taskNum])
        t_SCOOP_0 = time.time()
        react_species_task_results = list(map(rmg.reactionModel.react_speciesList, [newSpecies]*taskNum,
                                              corespeciesList_list))
        t_SCOOP_1 = time.time()
        t_SCOOP = t_SCOOP_1 - t_SCOOP_0

        ## processing the returned reactions by redirecting three main attributes:
        ## species, family, template
        t_getConst = []
        t_genRxn = []
        for idx in range(taskNum):
            reactions = react_species_task_results[idx][0]
            t_getConst.append(react_species_task_results[idx][1])
            t_genRxn.append(react_species_task_results[idx][2])
            for reaction in reactions:

                ### redirect family to family objects in root-worker
                ### print "Redirecting family and template for reactions..."
                reaction.family = families[reaction.family]

                ### redirect template to template objects in root-worker
                templateLabels = reaction.template
                redirect_template = []
                for label in templateLabels:
                    redirect_template.append(reaction.family.groups.entries[label])
                reaction.template = redirect_template

                ### de-IDize for species: convert ID into objects
                # print "De-IDizing for species in reactions..."
                reactants = []
                products = []
                pairs = []
                for reactant, product in reaction.pairs:
                    if isinstance(reactant, int):
                        reactant = rootSpeciesDict[reactant]
                    if isinstance(product, int):
                        product = rootSpeciesDict[product]
                    pairs.append((reactant, product))
                for reactant in reaction.reactants:
                    if isinstance(reactant, int):
                        reactant = rootSpeciesDict[reactant]
                    reactants.append(reactant)
                for product in reaction.products:
                    if isinstance(product, int):
                        product = rootSpeciesDict[product]
                    products.append(product)
                reaction.pairs = pairs
                reaction.products = products
                reaction.reactants = reactants

                ### redirect template for reaction.reverse
                # print "Redirecting family and template for reaction.reverse..."
                if hasattr(reaction, "reverse"):
                    reverseReaction = reaction.reverse
                    reverseReaction.family = families[reverseReaction.family]
                    reverseTemplateLabels = reverseReaction.template
                    redirect_reverseTemplate = []
                    for label in reverseTemplateLabels:
                        redirect_reverseTemplate.append(reverseReaction.family.groups.entries[label])
                    reverseReaction.template = redirect_reverseTemplate

            newReactions.extend(reactions)
        gc.collect()

        # output time profile for parallel mode
        entry = (caseNum, workerNum, t_setConst, t_SCOOP)
        for item in entry:
            entryStr += str(item)
            entryStr += '\t'
        time_profile_parallel = {"parallelMode": parallelMode, "t_setConst/sec": t_setConst, "t_SCOOP/sec": t_SCOOP,
                                 "t_getConst/sec": t_getConst, "t_genRxn/sec": t_genRxn, "max t_getConst/sec": max(t_getConst),
                                 "max t_genRxn/sec": max(t_genRxn)}
        with open("time_profile_parallel.json", 'a') as timeProfileFile:
            json.dump(time_profile_parallel, timeProfileFile)
            timeProfileFile.write('\n')

    # reaction generation done
    t2 = time.time()
    print "When parallelMode is set as {0}, the react time is {1}/s.".format(parallelMode, t2-t1)
    print "{0} reactions have been generated!".format(len(newReactions))

    # write parallel time profiling data
    if parallelMode:
        entryStr += str(t2-t1)
        entryStr += '\n'
    with open("time_profile_parallel.dat", 'a') as timeProfileFile:
            timeProfileFile.write(entryStr)
    # output time into times.json
    output_time = {"parallelMode": parallelMode, "time/sec": t2-t1, "rxns number":len(newReactions)}
    with open("times.json", 'a') as outputFile:
        json.dump(output_time, outputFile)
        outputFile.write('\n')


