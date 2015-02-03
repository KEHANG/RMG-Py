from scoop import futures, logger
from rmgpy.molecule.symmetryTest import *
import time
import json

def main():
    molecule = Molecule().fromSMILES('CC(C)CCC[C@@H](C)[C@H]1CC[C@@]2([H])[C@]3([H])CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C')#cholesterol
    # molecule = Molecule().fromSMILES('C(C(C1CCC1)(C2CCC2)C3CCC3)(C4CCC4)(C5CCC5)C6CCC6')#six-quares
    t_futures = 0
    t_ifelif = 0
    # parallel calculation
    n1 = time.time()
    tuple_parallel = calculateCyclicSymmetryNumber_parallel(molecule)
    n2 = time.time()
    t_parallel = (n2 - n1)*10**3
    symNum_par = tuple_parallel[0]
    t_eqGroup_par = tuple_parallel[1]
    t_eqBond_par = tuple_parallel[2]
    # serial calculation
    n3 = time.time()
    tuple_serial = calculateCyclicSymmetryNumber_serial(molecule)
    n4 = time.time()
    t_serial = (n4 - n3)*10**3
    symNum_ser = tuple_serial[0]
    t_eqGroup_ser = tuple_serial[1]
    t_eqBond_ser = tuple_serial[2]

    logger.info('CylicSymNum_parallel: {0}; CylicSymNum_serial: {1}\n'.format(symNum_par, symNum_ser))
    elapsedTime = {'t_parallel/ms': t_parallel, 't_eqGroup_par/ms': t_eqGroup_par, 't_eqBond_par/ms': t_eqBond_par,
                   't_serial/ms': t_serial, 't_eqGroup_ser/ms': t_eqGroup_ser, 't_eqBond_ser/ms': t_eqBond_ser}
    with open('times.json', 'a') as timeFile:
        json.dump(elapsedTime, timeFile)
        timeFile.write('\n')


if __name__ == "__main__":
    main()