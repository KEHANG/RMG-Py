from scoop import futures, logger
from rmgpy.molecule.symmetryTest import *
import time

def main():
    n1 = time.time()
    #print calculateAtomSymmetryNumber(molecule, atom)
    logger.info(testBondSymmetryNumberC8H18())
    n2 = time.time()
    t = (n2-n1)*10**3
    logger.info( 'Time: {0:.4f} milliseconds'.format(t))


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