from scoop import futures, logger
from rmgpy.molecule.symmetryTest import *
import time

def main():
    molecule = Molecule().fromSMILES('CC(C)(C)C(C)(C)C')
    t_parallel = 0
    t_non_parallel = 0
    t_futures = 0
    symmetryNumber = 1
    for atom1 in molecule.atoms:
        for atom2 in atom1.bonds:
            atomIdx1 = molecule.atoms.index(atom1)
            atomIdx2 = molecule.atoms.index(atom2)
            if atomIdx1 < atomIdx2:
                # parallel calculation
                n1 = time.time()
                tuple_s_t = calculateBondSymmetryNumber_parallel(molecule, atom1, atom2)
                n2 = time.time()
                symmetryNumber *= tuple_s_t[0]
                t_futures += tuple_s_t[1]
                t_parallel += (n2-n1)*10**3
                # for non-parallel calculation
                n1 = time.time()
                symmetryNumber *= calculateBondSymmetryNumber(molecule, atom1, atom2)
                n2 = time.time()
                t_non_parallel += (n2-n1)*10**3
    logger.info( 'Parallel Time: {0:.4f} milliseconds'.format(t_parallel))
    logger.info( 'Scoop Futures Time: {0:.4f} milliseconds'.format(t_futures))
    logger.info( 'Non_Parallel Time: {0:.4f} milliseconds'.format(t_non_parallel))

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