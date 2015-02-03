from scoop import futures, logger
from rmgpy.molecule.symmetryTest import *
import time
import json

def main():
    # molecule = Molecule().fromSMILES('CC(C)(C)C(C)(C)C')
    # molecule = Molecule().fromSMILES('CC(C)CCC[C@@H](C)[C@H]1CC[C@@]2([H])[C@]3([H])CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C')#cholesterol
    molecule = Molecule().fromSMILES('C(C(C1CCC1)(C2CCC2)C3CCC3)(C4CCC4)(C5CCC5)C6CCC6')#six-quares
    t_parallel = 0
    t_non_parallel = 0
    t_futures = 0
    t_ifelif = 0
    symmetryNumber1 = 1
    symmetryNumber2 = 1
    for atom1 in molecule.atoms:
        for atom2 in atom1.bonds:
            atomIdx1 = molecule.atoms.index(atom1)
            atomIdx2 = molecule.atoms.index(atom2)
            if atomIdx1 < atomIdx2:
                # parallel calculation
                n1 = time.time()
                # tuple_s_t = calculateBondSymmetryNumber_parallel(molecule, atom1, atom2)
                tuple_s_t1 = calculateBondSymmetryNumber_parallel(molecule, atom1, atom2)
                print 'tuple_s_t1: ', tuple_s_t1
                n2 = time.time()
                t_parallel += (n2-n1)*10**3
                symmetryNumber1 *= tuple_s_t1[0]
                t_futures += tuple_s_t1[1]

                # for non-parallel calculation
                n3 = time.time()
                tuple_s_t2 = calculateBondSymmetryNumber_ifelif(molecule, atom1, atom2)
                n4 = time.time()
                t_non_parallel += (n4-n3)*10**3
                symmetryNumber2 *= tuple_s_t2[0]
                t_ifelif += tuple_s_t2[1]
    # logger.info( 'Parallel Time: {0:.4f} milliseconds'.format(t_parallel))
    # logger.info( 'Scoop Futures Time: {0:.4f} milliseconds'.format(t_futures))
    # logger.info( 'Non_Parallel Time: {0:.4f} milliseconds'.format(t_non_parallel))
    logger.info('newSymNum: {0}; oldSymNum: {1}\n'.format(symmetryNumber1, symmetryNumber2))
    elapsedTime = {'t_parallel/ms': t_parallel, 't_futures/ms': t_futures,
                   't_non_parallel/ms': t_non_parallel, 't_ifelif/ms': t_ifelif}
    with open('times.json', 'a') as timeFile:
        json.dump(elapsedTime, timeFile)
        timeFile.write('\n')


def testBondSymmetryNumberEthane():
        """
        Test the Molecule.calculateBondSymmetryNumber() method.
        """
        molecule = Molecule().fromSMILES('CC')
        origSymmetryNumber = 1
        newSymmetryNumber = 1
        for atom1 in molecule.atoms:
            for atom2 in atom1.bonds:
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    origSymmetryNumber *= calculateBondSymmetryNumber(molecule, atom1, atom2)
                    newSymmetryNumber *= calculateBondSymmetryNumber_parallel(molecule, atom1, atom2)
        print 'origSymmetryNumber right?: ',origSymmetryNumber == 2
        print 'newSymmetryNumber right?: ', newSymmetryNumber ==2

def testBondSymmetryNumberC8H18():
        """
        Test the Molecule.calculateBondSymmetryNumber() method.
        """
        molecule = Molecule().fromSMILES('CC(C)(C)C(C)(C)C')
        origSymmetryNumber = 1
        newSymmetryNumber = 1
        for atom1 in molecule.atoms:
            for atom2 in atom1.bonds:
                if molecule.atoms.index(atom1) < molecule.atoms.index(atom2):
                    origSymmetryNumber *= calculateBondSymmetryNumber(molecule, atom1, atom2)
                    newSymmetryNumber *= calculateBondSymmetryNumber_parallel(molecule, atom1, atom2)
        print 'origSymmetryNumber right?: ',origSymmetryNumber == 2
        print 'newSymmetryNumber right?: ', newSymmetryNumber ==2

if __name__ == "__main__":
    main()