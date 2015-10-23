# global imports

import cython
import logging
import itertools
from collections import Counter

# local imports
try:
    import openbabel
except:
    pass
from rdkit import Chem

from rmgpy.molecule import element as elements

from .util import retrieveElementCount, VALENCES, ORDERS
from .inchi import AugmentedInChI, parse_H_layer, INCHI_PREFIX

from .pathfinder import \
 find_butadiene,\
 find_butadiene_end_with_charge,\
 find_allyl_end_with_charge

from .molecule import Atom, Bond, Molecule

# constants

BACKENDS = [
            'openbabel',
            'rdkit',
            ]

INCHI_LOOKUPS = {
            'H': '[H]',#RDkit was improperly handling the Hydrogen radical from InChI
            'He': '[He]',
        }
SMILES_LOOKUPS = {
            '[He]':# RDKit improperly handles helium and returns it in a triplet state
            """
            He
            multiplicity 1
            1 He u0 p1
            """
} 


def reset_lone_pairs_to_default(at):
    """Resets the atom's lone pair count to its default value."""

    bondorder = 0
    bonds = at.edges.values()
    for bond in bonds:
        bondorder += ORDERS[bond.order]
    
    at.lonePairs = (VALENCES[at.element.symbol] - bondorder - at.radicalElectrons - at.charge) / 2

def convert_unsaturated_bond_to_biradical(mol, inchi, u_indices):
    """
    Convert an unsaturated bond (double, triple) into a bond
    with a lower bond order (single, double), and give an unpaired electron
    to each of the neighboring atoms, with indices referring to the 1-based
    index in the InChI string.
    """
    cython.declare(u1=cython.int, u2=cython.int)
    cython.declare(atom1=Atom, atom2=Atom)
    cython.declare(b=Bond)

    combos = itertools.combinations(u_indices, 2)

    for u1, u2 in combos:
        atom1 = mol.atoms[u1 - 1] # convert to 0-based index for atoms in molecule
        atom2 = mol.atoms[u2 - 1] # convert to 0-based index for atoms in molecule
        if mol.hasBond(atom1, atom2):
            b = mol.getBond(atom1, atom2)
            if not b.isSingle():
                atom1.radicalElectrons += 1
                atom2.radicalElectrons += 1
                b.decrementOrder()

                u_indices.remove(u1)
                u_indices.remove(u2)

                return mol
            else:#maybe it's a mobile H-layer problem
                all_mobile_h_atoms_couples = parse_H_layer(inchi)
                """
                assume that O2=C1-O3H is the keto-enol system
                    1) find its partner (O2)
                    2) transfer H atom to partner (O2)
                    3) change bond order between partner and central carbon
                    4) add unpaired electrons to central carbon and original O.

                """
                if all_mobile_h_atoms_couples:
                    #find central atom:
                    central, original_atom, new_partner = find_mobile_h_system(mol, 
                        all_mobile_h_atoms_couples, [u1, u2])

                    # search Hydrogen atom and bond
                    hydrogen = None
                    for at, bond in original_atom.bonds.iteritems():
                        if at.number == 1:
                            hydrogen = at
                            mol.removeBond(bond)
                            break

                    new_h_bond = Bond(new_partner, hydrogen, order='S')
                    mol.addBond(new_h_bond)
                    
                    mol.getBond(central, new_partner).decrementOrder()

                    central.radicalElectrons += 1
                    original_atom.radicalElectrons += 1

                    u_indices.remove(u1)
                    u_indices.remove(u2)
                    return mol
        else:
            path = find_butadiene(atom1, atom2)
            if path is not None:
                atom1.radicalElectrons += 1
                atom2.radicalElectrons += 1
                # filter bonds from path and convert bond orders:
                bonds = path[1::2]#odd elements
                for bond in bonds[::2]:# even bonds
                    assert isinstance(bond, Bond)
                    bond.decrementOrder()
                for bond in bonds[1::2]:# odd bonds
                    assert isinstance(bond, Bond)
                    bond.incrementOrder()    

                u_indices.remove(u1)
                u_indices.remove(u2)

                return mol

    raise Exception('The indices {} did not refer to atoms that are connected in the molecule {}.'
        .format(u_indices, mol.toAdjacencyList()))    

def isUnsaturated(mol):
    """Does the molecule have a bond that's not single?
    
    (eg. a bond that is double or triple or beneze)"""
    cython.declare(atom1=Atom,
                   atom2=Atom,
                   bonds=dict,
                   bond=Bond)
    for atom1 in mol.atoms:
        bonds = mol.getBonds(atom1)
        for atom2, bond in bonds.iteritems():
            if not bond.isSingle():
                return True

    return False

def check_number_unpaired_electrons(mol):
    """Check if the number of unpaired electrons equals (m - 1)"""
    return mol.getNumberOfRadicalElectrons() == (mol.multiplicity - 1)        


def __fromSMILES(mol, smilesstr, backend):
    """Replace the Molecule `mol` with that given by the SMILES `smilesstr`
       using the backend `backend`"""
    if backend.lower() == 'rdkit':
        rdkitmol = Chem.MolFromSmiles(smilesstr)
        if rdkitmol is None:
            raise ValueError("Could not interpret the SMILES string {0!r}".format(smilesstr))
        fromRDKitMol(mol, rdkitmol)
        return mol
    elif backend.lower() == 'openbabel':
        parse_openbabel(mol, smilesstr, 'smi')
        return mol
    else:
        raise NotImplementedError('Unrecognized backend for SMILES parsing: {0}'.format(backend))

def __fromInChI(mol, inchistr, backend):
    """Replace the Molecule `mol` with that given by the InChI `inchistr`
       using the backend `backend`"""
    if backend.lower() == 'rdkit':
        rdkitmol = Chem.inchi.MolFromInchi(inchistr, removeHs=False)
        mol = fromRDKitMol(mol, rdkitmol)
        return mol 
    elif backend.lower() == 'openbabel':
        return parse_openbabel(mol, inchistr, 'inchi')
    else:
        raise NotImplementedError('Unrecognized backend for InChI parsing: {0}'.format(backend))


def __parse(mol, identifier, type_identifier, backend):
    """
    Parses the identifier based on the type of identifier (inchi/smi)
    and the backend used.
    
    First, look up the identifier in a dictionary to see if it can be processed
    this way.

    If not in the dictionary, parse it through the specified backed, 
    or try all backends.

    """

    if __lookup(mol, identifier, type_identifier) is not None:
        if isCorrectlyParsed(mol, identifier):
            return mol 

    for _backend in (BACKENDS if backend=='try-all' else [backend]):
        if type_identifier == 'smi':
            __fromSMILES(mol, identifier, _backend)
        elif type_identifier == 'inchi':
            __fromInChI(mol, identifier, _backend)
        else:
            raise NotImplementedError("Unknown identifier type {0}".format(type_identifier))

        if isCorrectlyParsed(mol, identifier):
            return mol
        else:
            logging.debug('Backend %s is not able to parse identifier %s', _backend, identifier)

    logging.error("Unable to correctly parse %s with backend %s", identifier, backend)
    raise Exception("Couldn't parse {0}".format(identifier))

def parse_openbabel(mol, identifier, type_identifier):
    """Converts the identifier to a Molecule using Openbabel."""
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(type_identifier, "smi")#SetInFormat(identifier) does not exist.
    obmol = openbabel.OBMol()
    obConversion.ReadString(obmol, identifier)
    obmol.AddHydrogens()
    obmol.AssignSpinMultiplicity(True)
    fromOBMol(mol, obmol)
    return mol


def isCorrectlyParsed(mol, identifier):
    """Check if molecule object has been correctly parsed."""
    conditions = []

    if mol.atoms:
        conditions.append(True)
    else:
        conditions.append(False)

    if 'InChI' in identifier:
        inchi_elementcount = retrieveElementCount(identifier)
        mol_elementcount = retrieveElementCount(mol)
        conditions.append(inchi_elementcount == mol_elementcount)

    return all(conditions)

def __lookup(mol, identifier, type_identifier):
    """
    Looks up the identifier and parses it the way we think is best.

    For troublesome inchis, we look up the smiles, and parse smiles.
    For troublesome smiles, we look up the adj list, and parse the adj list.

    """
    if type_identifier.lower() == 'inchi':
        try:
            smi = INCHI_LOOKUPS[identifier.split('/', 1)[1]]
            return mol.fromSMILES(smi)
        except KeyError:
            return None
    elif type_identifier.lower() == 'smi':
        try:
            adjList = SMILES_LOOKUPS[identifier]
            return mol.fromAdjacencyList(adjList)
        except KeyError:
            return None

def check(mol, aug_inchi) :
    """Check if molecule corresponds to the aug. inchi"""
    cython.declare(conditions=list,
                   inchi=str,
                   multi=cython.int,
                   )

    _, mult, __ = aug_inchi.inchi, aug_inchi.mult, aug_inchi.u_indices
    assert mult == mol.getRadicalCount() + 1,\
     'Multiplicity of molecule \n {0} does not correspond to aug. inchi {1}'.format(mol.toAdjacencyList(), aug_inchi)
    
    for at in mol.atoms:
        order = 0
        bonds = at.edges.values()
        for bond in bonds:
            order += ORDERS[bond.order]

        assert (order + at.radicalElectrons + 2*at.lonePairs + at.charge) == VALENCES[at.symbol],\
            'Valency for an atom of molecule \n {0} does not correspond to aug. inchi {1}'.format(mol.toAdjacencyList(), aug_inchi)

def correct_O_unsaturated_bond(mol, u_indices):
    """
    Searches for a radical or a charged oxygen atom connected to 
    a closed-shell carbon via an unsatured bond.

    Decrements the unsatured bond,
    transfers the unpaired electron from O to C or
    converts the charge from O to an unpaired electron on C, 
    increases the lone pair count of O to 2.

    Only do this once per molecule.
    """

    for at in mol.atoms:
        if at.isOxygen() and at.radicalElectrons == 1 and at.lonePairs == 1:
            bonds = mol.getBonds(at)
            oxygen = at
            for atom2, bond in bonds.iteritems():
                if bond.isTriple():
                    bond.decrementOrder()
                    oxygen.radicalElectrons -= 1
                    atom2.radicalElectrons += 1
                    oxygen.lonePairs += 1
                    return
        elif at.isOxygen() and at.charge == 1 and at.lonePairs == 1:
            bonds = mol.getBonds(at)
            oxygen = at

            start = oxygen
            # search for 3-atom-2-bond [X=X-X] paths
            paths = find_allyl_end_with_charge(start)
            for path in paths:    
                end = path[-1]
                start.charge += 1 if start.charge < 0 else -1
                end.charge += 1 if end.charge < 0 else -1
                start.lonePairs += 1
                # filter bonds from path and convert bond orders:
                bonds = path[1::2]#odd elements
                for bond in bonds[::2]:# even bonds
                    assert isinstance(bond, Bond)
                    bond.decrementOrder()
                for bond in bonds[1::2]:# odd bonds
                    assert isinstance(bond, Bond)
                    bond.incrementOrder()  
                return
            else:
                for atom2, bond in bonds.iteritems():
                    if not bond.isSingle() and atom2.charge == 0:
                        oxygen.charge -= 1
                        if (mol.atoms.index(atom2) + 1) in u_indices:
                            bond.decrementOrder()
                            atom2.radicalElectrons += 1
                            u_indices.remove(mol.atoms.index(atom2) + 1)
                        oxygen.lonePairs += 1
                        return

def fromInChI(mol, inchistr, backend='try-all'):
    """
    Convert an InChI string `inchistr` to a molecular structure. Uses 
    a user-specified backend for conversion, currently supporting
    rdkit (default) and openbabel.
    """

    mol.InChI = inchistr

    if INCHI_PREFIX in inchistr:
        return __parse(mol, inchistr, 'inchi', backend)
    else:
        return __parse(mol, INCHI_PREFIX + '/' + inchistr, 'inchi', backend)



def fromAugmentedInChI(mol, aug_inchi):
    """
    Creates a Molecule object from the augmented inchi.

    First, split off the multiplicity.
    Next, prepend the version layer to the inchi.
    Next, convert the inchi into a Molecule and
    set the multiplicity.

    Correct singlet one-center biradicals by replacing
    (u2, p0) by (u0, p1).

    Correct triplet two-center biradicals perceived as a
    zwitter ion by replacing (u0, c+/-1) by (u1, c0)

    Correct two-center triplet biradicals perceived as 
    a singlet double bond by replacing (u0)=(u0) by (u1)-(u1). 
    
    returns Molecule
    """

    if not isinstance(aug_inchi, AugmentedInChI):
        aug_inchi = AugmentedInChI(aug_inchi)

    mol = fromInChI(mol, aug_inchi.inchi)

    # multiplicity not specified in augmented InChI. Setting 
    if aug_inchi.mult == -1:
        logging.debug('Multiplicity not specified in augmented InChI.')
        logging.debug('Setting the multiplicity equal to the number of unpaired electrons + 1 of the parsed InChI.')
        mol.multiplicity = mol.getNumberOfRadicalElectrons() + 1
        return mol        

    mol.multiplicity = aug_inchi.mult

    #triplet to singlet conversion
    if mol.multiplicity == 1 and mol.getNumberOfRadicalElectrons() == 2:
        for at in mol.atoms:
            if at.radicalElectrons == 2:
                at.lonePairs = 1
                at.radicalElectrons = 0

    indices = aug_inchi.u_indices[:] if aug_inchi.u_indices is not None else []
    c = Counter(indices)
    for k,v in c.iteritems():
        atom = mol.atoms[k - 1]
        [indices.remove(k) for _ in range(atom.radicalElectrons)]


    contains_charge = sum([abs(at.charge) for at in mol.atoms]) != 0
    if mol.multiplicity >= 3 and not check_number_unpaired_electrons(mol) and contains_charge:
        fixCharge(mol, indices)

    # reset lone pairs                                
    for at in mol.atoms:
        reset_lone_pairs_to_default(at)

    # correct .O#C to O=C.
    correct_O_unsaturated_bond(mol, indices)


    # unsaturated bond to triplet conversion
    correct = check_number_unpaired_electrons(mol)

    unsaturated = isUnsaturated(mol)
    
    if not correct and not indices:
        raise Exception( 'Cannot correct {} based on {} by converting unsaturated bonds into unpaired electrons...'\
            .format(mol.toAdjacencyList(), aug_inchi))

    while not correct and unsaturated and len(indices) > 1:
        mol = convert_unsaturated_bond_to_biradical(mol, aug_inchi.inchi, indices)
        correct = check_number_unpaired_electrons(mol)
        unsaturated = isUnsaturated(mol)

    check(mol, aug_inchi)
    mol.updateAtomTypes()
    return mol

def fromSMILES(mol, smilesstr, backend='try-all'):
    """
    Convert a SMILES string `smilesstr` to a molecular structure. Uses 
    a user-specified backend for conversion, currently supporting
    rdkit (default) and openbabel.
    """
    return __parse(mol, smilesstr, 'smi', backend)


def fromSMARTS(mol, smartsstr):
    """
    Convert a SMARTS string `smartsstr` to a molecular structure. Uses
    `RDKit <http://rdkit.org/>`_ to perform the conversion.
    This Kekulizes everything, removing all aromatic atom types.
    """
    rdkitmol = Chem.MolFromSmarts(smartsstr)
    fromRDKitMol(mol, rdkitmol)
    return mol


def fromRDKitMol(mol, rdkitmol):
    """
    Convert a RDKit Mol object `rdkitmol` to a molecular structure. Uses
    `RDKit <http://rdkit.org/>`_ to perform the conversion.
    This Kekulizes everything, removing all aromatic atom types.
    """
    cython.declare(i=cython.int,
                   radicalElectrons=cython.int,
                   charge=cython.int,
                   lonePairs=cython.int,
                   number=cython.int,
                   order=cython.str,
                   atom=Atom,
                   atom1=Atom,
                   atom2=Atom,
                   bond=Bond)
    
    mol.vertices = []
    
    # Add hydrogen atoms to complete molecule if needed
    rdkitmol = Chem.AddHs(rdkitmol)
    Chem.rdmolops.Kekulize(rdkitmol, clearAromaticFlags=True)
    
    # iterate through atoms in rdkitmol
    for i in xrange(rdkitmol.GetNumAtoms()):
        rdkitatom = rdkitmol.GetAtomWithIdx(i)
        
        # Use atomic number as key for element
        number = rdkitatom.GetAtomicNum()
        element = elements.getElement(number)
            
        # Process charge
        charge = rdkitatom.GetFormalCharge()
        radicalElectrons = rdkitatom.GetNumRadicalElectrons()
        
        atom = Atom(element, radicalElectrons, charge, '', 0)
        mol.vertices.append(atom)
        
        # Add bonds by iterating again through atoms
        for j in xrange(0, i):
            rdkitatom2 = rdkitmol.GetAtomWithIdx(j + 1)
            rdkitbond = rdkitmol.GetBondBetweenAtoms(i, j)
            if rdkitbond is not None:
                order = ''
    
                # Process bond type
                rdbondtype = rdkitbond.GetBondType()
                if rdbondtype.name == 'SINGLE': order = 'S'
                elif rdbondtype.name == 'DOUBLE': order = 'D'
                elif rdbondtype.name == 'TRIPLE': order = 'T'
                elif rdbondtype.name == 'AROMATIC': order = 'B'
    
                bond = Bond(mol.vertices[i], mol.vertices[j], order)
                mol.addBond(bond)
    
    # Set atom types and connectivity values
    mol.update()
    mol.updateLonePairs()
    
    # Assume this is always true
    # There are cases where 2 radicalElectrons is a singlet, but
    # the triplet is often more stable, 
    mol.multiplicity = mol.getRadicalCount() + 1
    
    return mol

def fromOBMol(mol, obmol):
    """
    Convert a OpenBabel Mol object `obmol` to a molecular structure. Uses
    `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
    """
    # Below are the declared variables for cythonizing the module
    # cython.declare(i=cython.int)
    # cython.declare(radicalElectrons=cython.int, charge=cython.int, lonePairs=cython.int)
    # cython.declare(atom=Atom, atom1=Atom, atom2=Atom, bond=Bond)
    
    mol.vertices = []
    
    # Add hydrogen atoms to complete molecule if needed
    obmol.AddHydrogens()
    # TODO Chem.rdmolops.Kekulize(obmol, clearAromaticFlags=True)
    
    # iterate through atoms in obmol
    for obatom in openbabel.OBMolAtomIter(obmol):
        idx = obatom.GetIdx()#openbabel idx starts at 1!
        
        # Use atomic number as key for element
        number = obatom.GetAtomicNum()
        element = elements.getElement(number)
        # Process charge
        charge = obatom.GetFormalCharge()
        obatom_multiplicity = obatom.GetSpinMultiplicity()
        radicalElectrons =  obatom_multiplicity - 1 if obatom_multiplicity != 0 else 0
        
        atom = Atom(element, radicalElectrons, charge, '', 0)
        mol.vertices.append(atom)
    
    # iterate through bonds in obmol
    for obbond in openbabel.OBMolBondIter(obmol):
        order = 0
        # Process bond type
        oborder = obbond.GetBondOrder()
        if oborder == 1: order = 'S'
        elif oborder == 2: order = 'D'
        elif oborder == 3: order = 'T'
        elif obbond.IsAromatic() : order = 'B'

        bond = Bond(mol.vertices[obbond.GetBeginAtomIdx() - 1], mol.vertices[obbond.GetEndAtomIdx() - 1], order)#python array indices start at 0
        mol.addBond(bond)

    
    # Set atom types and connectivity values
    mol.updateConnectivityValues()
    mol.updateAtomTypes()
    mol.updateMultiplicity()
    mol.updateLonePairs()
    
    # Assume this is always true
    # There are cases where 2 radicalElectrons is a singlet, but
    # the triplet is often more stable, 
    mol.multiplicity = mol.getRadicalCount() + 1
    
    return mol

def fixCharge(mol, u_indices):
    """
    Fix molecules perceived as zwitterions that in reality are structures
    with multiple unpaired electrons.

    The simplest case converts atoms with a charge to atoms with one more
    unpaired electron.

    """
    if not u_indices:
        return

    # converting charges to unpaired electrons for atoms in the u-layer
    for at in mol.atoms:
        if at.charge != 0 and (mol.atoms.index(at) + 1) in u_indices:
            at.charge += 1 if at.charge < 0 else -1
            at.radicalElectrons += 1
            u_indices.remove(mol.atoms.index(at) + 1)

    # convert neighboring atoms (or delocalized paths) to unpaired electrons
    u_indices_copy = u_indices[:]
    for index in u_indices_copy:
        start = mol.atoms[index -1]

        # search for 4-atom-3-bond [X=X-X=X] paths
        path = find_butadiene_end_with_charge(start)
        if path is not None:    
            # we have found the atom we are looking for
            start.radicalElectrons += 1
            end = path[-1]
            end.charge += 1 if end.charge < 0 else -1
            end.lonePairs += 1
            # filter bonds from path and convert bond orders:
            bonds = path[1::2]#odd elements
            for bond in bonds[::2]:# even bonds
                assert isinstance(bond, Bond)
                bond.decrementOrder()
            for bond in bonds[1::2]:# odd bonds
                assert isinstance(bond, Bond)
                bond.incrementOrder()  
            u_indices.remove(mol.atoms.index(start) + 1)
            continue

        # search for 3-atom-2-bond [X=X-X] paths
        paths = find_allyl_end_with_charge(start)
        from rmgpy.data.kinetics.family import ReactionRecipe

        for path in paths:
            # label atoms so that we can use the labels in the actions of the recipe
            for i, at in enumerate(path[::2]):
                assert isinstance(at, Atom)
                at.label = str(i)
            # we have found the atom we are looking for
            fix_charge_recipe = ReactionRecipe()
            fix_charge_recipe.addAction(['GAIN_RADICAL', start.label, 1])

            end = path[-1]
            end_original_charge = end.charge
          
            # filter bonds from path and convert bond orders:
            bonds = path[1::2]#odd elements
            for bond in bonds[::2]:# even bonds
                assert isinstance(bond, Bond)
                fix_charge_recipe.addAction(['CHANGE_BOND', bond.atom1.label, -1, bond.atom2.label])
            for bond in bonds[1::2]:# odd bonds
                assert isinstance(bond, Bond)
                fix_charge_recipe.addAction(['CHANGE_BOND', bond.atom1.label, 1, bond.atom2.label])

            end.charge += 1 if end.charge < 0 else -1
            fix_charge_recipe.applyForward(mol, update=False)

            if check_bond_order_oxygen(mol):
                u_indices.remove(mol.atoms.index(start) + 1)
                # unlabel atoms so that they never cause trouble downstream
                for i, at in enumerate(path[::2]):
                    assert isinstance(at, Atom)
                    at.label = ''
                break
            else:
                fix_charge_recipe.applyReverse(mol, update=False)
                end.charge = end_original_charge

                # unlabel atoms so that they never cause trouble downstream
                for i, at in enumerate(path[::2]):
                    assert isinstance(at, Atom)
                    at.label = ''

                continue # to next path

            
        continue # to next index in u-layer

    # fix adjacent charges
    for at in mol.atoms:
        if at.charge != 0:
            for neigh, bond in at.bonds.iteritems():
                if neigh.charge != 0:
                    bond.incrementOrder()
                    at.charge += 1 if at.charge < 0 else -1
                    neigh.charge += 1 if neigh.charge < 0 else -1

def find_mobile_h_system(mol, all_mobile_h_atoms_couples, test_indices):
    dummy = test_indices[:]

    for mobile_h_atom_couple in all_mobile_h_atoms_couples:
        for test_index in test_indices:
            if test_index in mobile_h_atom_couple:
                original_atom = test_index
                dummy.remove(test_index)
                mobile_h_atom_couple.remove(test_index)
                new_partner = mobile_h_atom_couple[0]
                central = dummy[0]
                return mol.atoms[central - 1], mol.atoms[original_atom - 1], mol.atoms[new_partner - 1]

    raise Exception('We should always have found the mobile-H system. All mobile H couples: {}, test indices: {}'
        .format(all_mobile_h_atoms_couples, test_indices))

def check_bond_order_oxygen(mol):
    """Check if total bond order of oxygen atoms is smaller than 4."""
    from rmgpy.molecule.util import ORDERS

    for at in mol.atoms:
        if at.number == 8:
            order = sum([ORDERS[b.order] for _, b in at.bonds.iteritems()])
            not_correct = order >= 4
            if not_correct:
                return False

    return True
    
