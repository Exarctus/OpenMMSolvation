from simtk.openmm.app import *
from simtk.openmm import *
from simtk import unit
from simtk.openmm import app
from simtk.openmm.app import PDBFile, Modeller, PDBFile
from mdtraj.reporters import NetCDFReporter   
from openmmtools import alchemy
import numpy as np
from time import time

from openforcefield.topology import Molecule, Topology
from openmmforcefields.generators import SystemGenerator
from openmmtools.forces import find_forces
from openmmtools.constants import ONE_4PI_EPS0
import copy
from openmmtools import forcefactories
from openmmtools import forces


def remake_system(system, solute_indexes, cutoff=9.0 * unit.angstroms, switching_distance=8.0 * unit.angstroms):
    
    new_system = copy.deepcopy(system)
    
    alchemical_atoms = set(solute_indexes)
    chemical_atoms = set(range(system.getNumParticles())) - alchemical_atoms
    
    force_idx, reference_force = find_forces(new_system, openmm.NonbondedForce, only_one=True)
    
    nonbonded_force = copy.deepcopy(reference_force)
    
    epsilon_solvent = nonbonded_force.getReactionFieldDielectric()
    
    k_rf = cutoff ** (-3) * ((epsilon_solvent - 1) / (2 * epsilon_solvent + 1))
    k_rf = k_rf.value_in_unit_system(unit.md_unit_system)  
    
    electrostatics_function = 'ONE_4PI_EPS0*chargeprod*(r^(-1) +k_rf*r^2);'
    electrostatics_function += 'k_rf = {k_rf};'.format(k_rf=k_rf)
    electrostatics_function += 'ONE_4PI_EPS0 = %.16e;' % (ONE_4PI_EPS0)
    electrostatics_function += 'chargeprod = charge1*charge2'
    electrostatics = CustomNonbondedForce(electrostatics_function)
    electrostatics.addPerParticleParameter('charge') 

    electrostatics_bond_function = 'ONE_4PI_EPS0*chargeprod*(r^(-1) +k_rf*r^2);'
    electrostatics_bond_function += 'k_rf = {k_rf};'.format(k_rf=k_rf)
    electrostatics_bond_function += 'ONE_4PI_EPS0 = %.16e;' % (ONE_4PI_EPS0)
    electrostatics_bond = CustomBondForce(electrostatics_bond_function)
    electrostatics_bond.addPerBondParameter("chargeprod") 
                
    lj_function = '4.0*epsilon*x*(x-1.0); x = (sigma/r)^6;'
    lj_function += 'sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)'
    lj = CustomNonbondedForce(lj_function)
    lj.addPerParticleParameter('sigma')
    lj.addPerParticleParameter('epsilon')
    
    lj_bond_function = '4.0*epsilon*x*(x-1.0); x = (sigma/r)^6;'
    sterics_bond = CustomBondForce(lj_bond_function)
    
    sterics_bond.addPerBondParameter("sigma")
    sterics_bond.addPerBondParameter("epsilon")
    
    for ind in range(nonbonded_force.getNumParticles()):

        [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(ind)

        if sigma / unit.nanometer == 0.0:
            newsigma = 0.3 * unit.nanometer 
        else:
            newsigma = sigma

        lj.addParticle([sigma, epsilon])
        electrostatics.addParticle([charge])
    
    for particle_index in range(nonbonded_force.getNumParticles()):
        [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(particle_index)

        if particle_index in alchemical_atoms:
            # not necessary as this force gets deleted, but can't hurt.
            nonbonded_force.setParticleParameters(particle_index, abs(0.0 * charge), sigma, abs(0 * epsilon))
                    
    for ind in range(nonbonded_force.getNumExceptions()):
        [iatom, jatom, chargeprod, sigma, epsilon] = nonbonded_force.getExceptionParameters(ind)
        
        lj.addExclusion(iatom, jatom)
        electrostatics.addExclusion(iatom, jatom)
        
        is_exception_epsilon = abs(epsilon.value_in_unit_system(unit.md_unit_system)) > 0.0
        is_exception_chargeprod = abs(chargeprod.value_in_unit_system(unit.md_unit_system)) > 0.0
        
        both_alchemical = iatom in alchemical_atoms and jatom in alchemical_atoms
        at_least_one_alchemical = iatom in alchemical_atoms or jatom in alchemical_atoms
        only_one_alchemical = at_least_one_alchemical and not both_alchemical
        
        if both_alchemical:
            if is_exception_epsilon:
                sterics_bond.addBond(iatom, jatom, [sigma, epsilon])
            if is_exception_chargeprod:
                electrostatics_bond.addBond(iatom, jatom, [chargeprod])
                
        elif only_one_alchemical:
            if is_exception_epsilon:
                sterics_bond.addBond(iatom, jatom, [sigma, epsilon])
            if is_exception_chargeprod:
                electrostatics_bond.addBond(iatom, jatom, [chargeprod])

        if at_least_one_alchemical:
            # not necessary as this force gets deleted, but can't hurt.
            nonbonded_force.setExceptionParameters(ind, iatom, jatom, abs(0.0 * chargeprod), sigma, abs(0.0 * epsilon))
        
    electrostatics.addInteractionGroup(alchemical_atoms, alchemical_atoms)
    electrostatics.addInteractionGroup(alchemical_atoms, chemical_atoms)
    electrostatics.addInteractionGroup(chemical_atoms, chemical_atoms)
    
    lj.addInteractionGroup(alchemical_atoms, alchemical_atoms)
    lj.addInteractionGroup(alchemical_atoms, chemical_atoms)
    lj.addInteractionGroup(chemical_atoms, chemical_atoms)
    
    electrostatics.setCutoffDistance(cutoff)
    electrostatics.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    electrostatics.setUseLongRangeCorrection(False)
    
    if (switching_distance is not None):
        electrostatics.setUseSwitchingFunction(True)
        electrostatics.setSwitchingDistance(switching_distance)
    else: 
        electrostatics.setUseSwitchingFunction(False)
        
    electrostatics.setForceGroup(0)
    electrostatics_bond.setForceGroup(0)
        
    new_system.removeForce(force_idx)
    
    new_system.addForce(electrostatics)
    new_system.addForce(electrostatics_bond)
    
    lj.setCutoffDistance(cutoff)
    lj.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    lj.setUseLongRangeCorrection(False) 
    
    if (switching_distance is not None):
        lj.setUseSwitchingFunction(True)
        lj.setSwitchingDistance(switching_distance)
    else: 
        lj.setUseSwitchingFunction(False)
        
    lj.setForceGroup(0)
    sterics_bond.setForceGroup(0)
    
    new_system.addForce(lj)
    new_system.addForce(sterics_bond)
    
    return new_system


def collect_solute_indexes(topology):
    soluteIndices = []
    for res in topology.residues():
        resname = res.name.upper()
        if (resname != 'HOH' and resname != 'WAT'and resname != 'CL'and resname != 'NA'):
            for atom in res.atoms():
                soluteIndices.append(atom.index)
    return soluteIndices


platform = openmm.Platform.getPlatformByName('OpenCL')
platform.setPropertyDefaultValue('Precision', 'mixed')

'''
---SYSTEM PREPARATION---
    setup AM1-BCC charges for the solute, add solvent, set non-bonded method etc
'''
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

new_sys = remake_system(system, collect_solute_indexes(modeller.topology), cutoff=9.0 * unit.angstroms, switching_distance=7.5 * unit.angstroms)

'''
---FINISHED SYSTEM PREPARATION---
'''
    
# new_sys.addForce(MonteCarloBarostat(1 * unit.bar, 298.15 * unit.kelvin))
integrator = LangevinIntegrator(298.15 * unit.kelvin, 1.0 / unit.picoseconds, 0.002 * unit.picoseconds)
# integrator.setIntegrationForceGroups(set([0]))

simulation = app.Simulation(modeller.topology, new_sys, integrator, platform)
simulation.context.setPositions(modeller.positions)

simulation.minimizeEnergy()

state = simulation.context.getState(getPositions=True)
PDBFile.writeFile(modeller.topology, state.getPositions(), file=open("equil.pdb", "w"))

simulation.reporters.append(StateDataReporter('data.txt', 100, step=True, potentialEnergy=True, temperature=True, density=True , volume=True))
simulation.reporters.append(NetCDFReporter('output.nc', 100))

start = time()
simulation.step(10000)
end = time()
print (end - start)

