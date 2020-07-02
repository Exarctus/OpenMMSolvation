# OpenMM CustomNonbondedForce/Alchemy.py Tests


This repository was originally intended to calculate the free energy of solvation of ethanol in water at standard temperature and pressure using openmm, openmmtools and openmmforcefields. First, however, it seems there's an issue with how CustomNonbondedForce implements/uses switching functions to remove energy discontinuity at the cutoff. 

##Requirements
openmm
openmmtools
openmmforcefields

## Using Alchemy.py and Implementational Details

To use the code, simply type:

```python
python3 alchemy_test.py
```
This uses alchemy.py as the driver for decoupling solvent and solute interactions. 

The solvated ethanol-water system is first prepared using SystemGenerator. A 12A TIP4PEW waterbox with nonbonded interactions cutoff set to 9A was created. Long range electrostatic interactions are evaluated with the Reaction Field method using a switching distance of 7.5 angstrom (seems to be common in literature).

```python
ligand_mol = Molecule.from_file('ethanol.sdf', file_format='sdf')

forcefield_kwargs = {'constraints': app.HBonds, 'rigidWater': True, 'removeCMMotion': True, 'hydrogenMass': 4 * unit.amu }

system_generator = SystemGenerator(
    forcefields=['amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml'],
    small_molecule_forcefield='gaff-2.11',
    molecules=[ligand_mol],
    forcefield_kwargs=forcefield_kwargs)

ligand_pdb = PDBFile('ethanol.pdb')

modeller = Modeller(ligand_pdb.topology, ligand_pdb.positions)

modeller.addSolvent(system_generator.forcefield, model='tip3p', padding=12.0 * unit.angstroms)

system = system_generator.forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffPeriodic,
        nonbondedCutoff=9.0 * unit.angstroms, constraints=HBonds, switchDistance=7.5 * unit.angstroms)
```

The alchemical system is then prepared using openmmtools' alchemy.py code:

```python
def collect_solute_indexes(topology):
    soluteIndices = []
    for res in topology.residues():
        resname = res.name.upper()
        if (resname != 'HOH' and resname != 'WAT'and resname != 'CL'and resname != 'NA'):
            for atom in res.atoms():
                soluteIndices.append(atom.index)
    return soluteIndices
    
factory = AbsoluteAlchemicalFactory(consistent_exceptions=False)
alchemical_region = AlchemicalRegion(alchemical_atoms=collect_solute_indexes(modeller.topology), annihilate_sterics=False, annihilate_electrostatics=False)
alchemical_system = factory.create_alchemical_system(system, alchemical_region)

alchemical_state = AlchemicalState.from_system(alchemical_system)
alchemical_state.lambda_electrostatics = 1.0
alchemical_state.lambda_sterics = 1.0
alchemical_state.apply_to_system(alchemical_system)
```

Here, annihiliate_electrostatics/sterics are set to false as we ideally only want to decouple the interactions. The system is then run in NVT for a short period of time **at the fully interacting end point**.

```python
# new_sys.addForce(MonteCarloBarostat(1 * unit.bar, 298.15 * unit.kelvin))
integrator = LangevinIntegrator(298.15 * unit.kelvin, 1.0 / unit.picoseconds, 0.002 * unit.picoseconds)
# integrator.setIntegrationForceGroups(set([0]))

simulation = app.Simulation(modeller.topology, alchemical_system, integrator, platform)
simulation.context.setPositions(modeller.positions)

simulation.minimizeEnergy()

state = simulation.context.getState(getPositions=True)
PDBFile.writeFile(modeller.topology, state.getPositions(), file=open("equil.pdb", "w"))

simulation.reporters.append(StateDataReporter('data.txt', 100, step=True, potentialEnergy=True, temperature=True, density=True , volume=True))
simulation.reporters.append(NetCDFReporter('output.nc', 100))

```

Opening the resulting MD files in VMD (vmd equil.pdb output.nc) shows that the atoms are more or less frozen during the course of dynamics. Here, the system should be completely coupled, however there seems something weird going on, so I tried to remake my own version of alchemy.py.


## Using custom_force_test.py and Implementational Details

To use the code, simply type:

```python
python3 custom_force_test.py
```

this code is nearly identical to the above, however rather than using alchemy.py it includes a function which redefines the solute-solvent, as well as the solvent-solvent and solute-solute interactions:

```python
new_sys = remake_system(system, collect_solute_indexes(modeller.topology), cutoff=9.0 * unit.angstroms, switching_distance=7.5 * unit.angstroms)
```

this function first defines reaction-field electrostatics and lennard-jones interactions, adds the corresponding exclusions to each, with exceptions being modelled by CustomBondForces (not sure if there's a better way to do this?). Interaction groups are then used to define solute- and solvent- specific interaction sets, and finally the cutoff, switching and nonbonded method types are set for each.

Finally, the original NonbondedForce gets deleted, and replaced by all the custom forces representing the missing non-bonded interactions.

Using a switching function starting at 7.5 angstroms leads to identical issues as with alchemy.py.

**It should be noted, that running dynamics with a normal NonbondedForce (not custom) seems to be perfectly capable of handling SwitchingFunctions, it seems it's only when using CustomNonbondedForce that the issues arise**. 



