# global imports

cimport element as elements

# no .pxd files for these:
#from .util cimport retrieveElementCount, VALENCES, ORDERS
#from .inchi cimport AugmentedInChI, compose_aug_inchi_key, compose_aug_inchi, INCHI_PREFIX, MULT_PREFIX, U_LAYER_PREFIX

from .molecule cimport Atom, Bond, Molecule

cpdef list BACKENDS
cpdef dict INCHI_LOOKUPS
cpdef dict SMILES_LOOKUPS
cpdef dict _known_smiles_molecules
cpdef _known_smiles_radicals

cpdef reset_lone_pairs_to_default(Atom at)

cdef Molecule convert_unsaturated_bond_to_biradical(Molecule mol, list u_indices)

cpdef bint isUnsaturated(Molecule mol)
    
cpdef bint check_number_unpaired_electrons(Molecule mol)

cdef Molecule __fromSMILES(Molecule mol, str smilesstr, str backend)

cdef Molecule __fromInChI(Molecule mol, str inchistr, str backend)

cdef __parse(Molecule mol, str identifier, str type_identifier, str backend)

cpdef parse_openbabel(Molecule mol, str identifier, str type_identifier)

cpdef isCorrectlyParsed(Molecule mol, str identifier)

cdef __lookup(Molecule mol, str identifier, str type_identifier)

cpdef isZwitterIon(Molecule mol)
   
cpdef check(Molecule mol, aug_inchi)

cpdef correct_O_triple_bond(Molecule mol)

cpdef fromInChI(Molecule mol, str inchistr, backend=*)

cpdef Molecule fromAugmentedInChI(Molecule mol, aug_inchi)

cpdef fromSMILES(Molecule mol, str smilesstr, str backend=*)


cpdef fromSMARTS(Molecule mol, str smartsstr)
    
cpdef Molecule fromRDKitMol(Molecule mol, object rdkitmol)

cpdef fromOBMol(Molecule mol, object obmol)

cpdef str toSMARTS(Molecule mol)

cpdef str toSMILES(Molecule mol)

cpdef toOBMol(Molecule mol)

cpdef toRDKitMol(Molecule mol, bint removeHs=*, bint returnMapping=*, bint sanitize=*)

cpdef str toInChI(Molecule mol)

cpdef createULayer(Molecule mol)

# returns an AugmentedInChI but there's no pxd file for that yet
cpdef toAugmentedInChI(Molecule mol)

cpdef str toInChIKey(Molecule mol)

cpdef str toAugmentedInChIKey(Molecule mol)

cpdef str createMultiplicityLayer(int multiplicity)

cpdef fixZwitter(Molecule mol)

cpdef sortAtoms(Molecule mol)

cpdef moveHs(Molecule mol)

cpdef updateAtomConnectivityValues(Molecule mol)

cpdef Molecule normalize(Molecule mol)

cpdef list get_unpaired_electrons(Molecule mol)

cpdef list find_delocalized_path(Atom start, Atom end)

cpdef list findAllylPaths(list existing_path)
