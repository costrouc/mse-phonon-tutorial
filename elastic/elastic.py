import json

from pymatgen import Structure, Lattice
from pymatgen.analysis.elasticity.strain import DeformedStructureSet, Strain
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.analysis.elasticity.elastic import ElasticTensor

def create_structure():
    # Create MgO Conventional Unit Cell
    a = 4.292434060
    lattice = Lattice.from_parameters(a, a, a, 90, 90, 90)
    symbols = ['Mg', 'O']
    positions = [[0, 0, 0], [0.5, 0.5, 0.5]]
    mgo_conventional = Structure.from_spacegroup(225, lattice, symbols, positions)
    return mgo_conventional

def write_espresso_file(structure, index):
    input_template = """
&control
    calculation = 'relax'
    pseudo_dir = '../pseudo/'
    outdir = './output-{index}/'
    forc_conv_thr = 1.d-7
    tprnfor = .true.
    tstress = .true.
 /
 &system
    ibrav = 0
    nat = 8
    ntyp = 2
    ecutwfc = 70.0
 /
 &electrons
    diagonalization='david'
    conv_thr = 1.d-9
 /
 &ions
    ion_dynamics = 'bfgs'
 /
ATOMIC_SPECIES
 Mg   24.30500   Mg.pbe-nsp-bpaw.UPF
  O   15.99940   O.pbe-kjpaw.UPF
K_POINTS automatic
 6 6 6 1 1 1""".format(index=index)
    cell_parameters = [' {:0.8f} {:0.8f} {:0.8f}'.format(*vector) for vector in structure.lattice.matrix]
    atom_positions = [' {} {:0.8f} {:0.8f} {:0.8f}'.format(atom.specie.name, *atom.frac_coords) for atom in structure]
    input_text = [
        input_template, 
        'CELL_PARAMETERS angstrom',
        '\n'.join(cell_parameters),
        'ATOMIC_POSITIONS crystal',
        '\n'.join(atom_positions),
        '\n'
    ]
    filename = 'relax-%d.in' % index
    with open(filename, 'w') as f:
        f.write('\n'.join(input_text))

def create_deformed_structures(structure):
    # Create Deformed Structures
    max_normal_deformation = 0.01
    max_shear_deformation = 0.02 # Default value of 0.08 is TOO high
    num_normal = 4
    num_shear = 4

    deformation_set = DeformedStructureSet(structure, 
                                           nd=max_normal_deformation,
                                           ns=max_shear_deformation,
                                           num_norm=num_normal,
                                           num_shear=num_shear)

    data = []
    for i, deformation in enumerate(deformation_set.deformations):
        deformed_structure = deformation.apply_to_structure(structure)
        strain = Strain.from_deformation(deformation).tolist()
        data.append({'index': i+1, 'structure': deformed_structure , 'strain': strain})
    return data

def read_stresses_from_output(filename):
    with open(filename) as f:
        data = f.read()

    stress = []
    for line in data[data.rfind('total   stress'):].split('\n')[1:4]:
        stress.append(list(map(float, line.split()[3:])))
    return Stress(stress) / -10 # unit conversion

def print_elastic_analysis(elastic):
    print('Stiffness Tensor')
    for row in elastic.voigt:
        print('{:+8.1f} {:+8.1f} {:+8.1f} {:+8.1f} {:+8.1f} {:+8.1f}\n'.format(*row))
    print('\n')
    print('Shear Modulus G_V', elastic.g_voigt)
    print('Shear Modulus G_R', elastic.g_reuss)
    print('Shear Modulus G_vrh', elastic.g_vrh)

    print('Bulk Modulus K_V', elastic.k_voigt)
    print('Bulk Modulus K_R', elastic.k_reuss)
    print('Bulk Modulus K_vrh', elastic.k_vrh)

    print('Elastic Anisotropy', elastic.universal_anisotropy)
    print('Poisons Ration', elastic.homogeneous_poisson)

    
def build():
    structure = create_structure()
    data = create_deformed_structures(structure)
    for deformed_structure in data:
        # Hackish way to write input file
        write_espresso_file(deformed_structure['structure'],
                            deformed_structure['index']) 

def analyze():
    structure = create_structure()
    data = create_deformed_structures(structure)
    for deformed_structure in data:
        filename = 'relax-%d.out' % deformed_structure['index']
        stress = read_stresses_from_output(filename)
        deformed_structure['stress'] = stress

    stresses = [s['stress'] for s in data]
    strains = [s['strain'] for s in data]
    elastic = ElasticTensor.from_strain_stress_list(strains, stresses)
    print_elastic_analysis(elastic)


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2 or sys.argv[1] not in ['build', 'analysis']:
        print('must be run with either build or analyze')
        exit()

    option = sys.argv[1]
    if option == 'build':
        build()
    elif option == 'analysis':
        analyze()
