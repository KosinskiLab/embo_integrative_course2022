#!/usr/bin/env python3

from __future__ import print_function

# Imports needed to use ProtocolOutput
import IMP.pmi.mmcif
import ihm

import IMP
import IMP.core
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.tools
import IMP.pmi.samplers
import RMF
import IMP.rmf

import IMP.pmi.macros
import IMP.pmi.topology
import tutorial_util

import sys

import warnings
warnings.filterwarnings('ignore')


# Then setup the relevant paths of the input files

# In[2]:
try:
    import IMP.mpi
    print('ReplicaExchange: MPI was found. Using Parallel Replica Exchange')
    rex_obj = IMP.mpi.ReplicaExchange()
except ImportError:
    print('ReplicaExchange: Could not find MPI. Using Serial Replica Exchange')
    rex_obj = IMP.pmi.samplers._SerialReplicaExchange()

replica_number = rex_obj.get_my_index()


output_prefix=""
step=int(1)

datadirectory = "../data/"
topology_file = datadirectory+"topology_poliii.cryoem.txt"
target_gmm_file=datadirectory+'%d_imp.gmm'%(step)
output_dir='.'
output_directory='%s/%s_output_%d/'%(output_dir,output_prefix, step)



# Initialize IMP model
m = IMP.Model()

# Read in the topology file.  
# Specify the directory wheere the PDB files, fasta files and GMM files are
topology = IMP.pmi.topology.TopologyReader(topology_file,
                                  pdb_dir=datadirectory,
                                  fasta_dir=datadirectory,
                                  gmm_dir=datadirectory)




# Use the BuildSystem macro to build states from the topology file
bs = IMP.pmi.macros.BuildSystem(m)

# Record the modeling protocol to an mmCIF file
po = IMP.pmi.mmcif.ProtocolOutput(None)
bs.system.add_protocol_output(po)
po.system.title = "Modeling of RNA Pol III"
# Add publication
po.system.citations.append(ihm.Citation.from_pubmed_id(25161197))

# Each state can be specified by a topology file.
bs.add_state(topology)




root_hier, dof = bs.execute_macro(max_rb_trans=4.0,
                                  max_rb_rot=0.3,
                                  max_bead_trans=4.0,
                                  max_srb_trans=4.0,
                                  max_srb_rot=0.3)





import IMP.pmi.plotting
import IMP.pmi.plotting.topology

IMP.pmi.plotting.topology.draw_component_composition(dof)

fixed_particles=[]
for prot in ["ABC23"]:
    fixed_particles+=IMP.atom.Selection(root_hier,molecule=prot).get_selected_particles()


# Fix the Corresponding Rigid movers and Super Rigid Body movers using dof
# The flexible beads will still be flexible (fixed_beads is an empty list)!
fixed_beads,fixed_rbs=dof.disable_movers(fixed_particles,
                                         [IMP.core.RigidBodyMover,IMP.pmi.TransformMover])




if step==1:
    # Shuffle the rigid body configuration of only the molecules we are interested in (Rpb4 and Rpb7)
    # but all flexible beads will also be shuffled.
    IMP.pmi.tools.shuffle_configuration(root_hier,
                                        max_translation=300,
                                        verbose=True,
                                        cutoff=5.0,
                                        niterations=100)
                                        #excluded_rigid_bodies=fixed_rbs,


else:
    rh_ref = RMF.open_rmf_file_read_only('seed_%d.rmf3'%(step-1))
    IMP.rmf.link_hierarchies(rh_ref, [root_hier])
    IMP.rmf.load_frame(rh_ref, RMF.FrameID(replica_number))

outputobjects = [] # reporter objects...output is included in the stat file



# Connectivity keeps things connected along the backbone (ignores if inside same rigid body)
mols = IMP.pmi.tools.get_molecules(root_hier)
for mol in mols:
    molname=mol.get_name()
    IMP.pmi.tools.display_bonds(mol)
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol,scale=2.0)
    cr.add_to_model()
    cr.set_label(molname)
    outputobjects.append(cr)


# #### Excluded Volume Restraint <a name="Excluded_Volume_Restraint_4"></a>

# In[12]:


ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                         included_objects=root_hier,
                                         resolution=10)
ev.add_to_model()         # add to scoring function
outputobjects.append(ev)  # add to output


# #### Crosslinks - dataset 1 <a name="Crosslink_1_4"></a>
# 


# We then initialize a CrossLinkDataBase that uses a keywords converter to map column to information.
# The required fields are the protein and residue number for each side of the crosslink.
xldbkwc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkwc.set_protein1_key("Protein1")
xldbkwc.set_protein2_key("Protein2")
xldbkwc.set_residue1_key("AbsPos1")
xldbkwc.set_residue2_key("AbsPos2")
xldbkwc.set_id_score_key("ld-Score")

xl1 = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
xl1.create_set_from_file(datadirectory+'FerberKosinski2016_apo.csv')
xl1.set_name("APO")

xl2 = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
xl2.create_set_from_file(datadirectory+'FerberKosinski2016_DNA.csv')
xl2.set_name("DNA")

# Append the xl2 dataset to the xl1 dataset to create a larger dataset
xl1.append_database(xl2)

# Rename one protien name
xl1.rename_proteins({"ABC14.5":"ABC14_5"})

# Create 3 confidence classes
xl1.classify_crosslinks_by_score(3)



# Now, we set up the restraint.
xl1rest = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                                   root_hier=root_hier,  # The root hierarchy
                                   CrossLinkDataBase=xl1,# The XLDB defined above
                                   length=21.0,          # Length of the linker in angstroms
                                   slope=0.002,           # A linear term that biases XLed
                                                         # residues together
                                   resolution=1.0,       # Resolution at which to apply the restraint. 
                                                         # Either 1 (residue) or 0 (atomic)
                                   label="XL",         # Used to label output in the stat file
                                   weight=10.)            # Weight applied to all crosslinks 
                                                         # in this dataset
xl1rest.add_to_model()
outputobjects.append(xl1rest)

# Electron Microscopy Restraint
#  The GaussianEMRestraint uses a density overlap function to compare model to data
#   First the EM map is approximated with a Gaussian Mixture Model (done separately)
#   Second, the components of the model are represented with Gaussians (forming the model GMM)

import IMP.bayesianem
import IMP.bayesianem.restraint
# First, get the model density objects that will be fitted to the EM density.
densities = IMP.atom.Selection(root_hier, representation_type=IMP.atom.DENSITIES).get_selected_particles()
gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(densities,
                                                 target_fn=target_gmm_file,
                                                 scale_target_to_mass=True,
                                                 slope=0.01,
                                                 target_radii_scale=3.0,
                                                 target_is_rigid_body=False)

gem.add_to_model()
gem.set_label("Total")
outputobjects.append(gem)

# total number of saved frames
num_frames = 2

# This object defines all components to be sampled as well as the sampling protocol
mc1=IMP.pmi.macros.ReplicaExchange0(m,
              root_hier=root_hier,                         # The root hierarchy
              monte_carlo_sample_objects=dof.get_movers()+xl1rest.get_movers(), # All moving particles and parameters
              rmf_output_objects=outputobjects,            # Objects to put into the rmf file
              crosslink_restraints=[xl1rest],      # allows XLs to be drawn in the RMF files
              monte_carlo_temperature=1.0,
              replica_exchange_minimum_temperature=1.0,
              replica_exchange_maximum_temperature=2.5,
              simulated_annealing=False,
              number_of_best_scoring_models=0,
              monte_carlo_steps=10,
              number_of_frames=num_frames,
              #save_coordinates_mode="25th_score",
              global_output_directory=output_directory,
              test_mode=True)

# Start Sampling
mc1.execute_macro()

po.finalize()

s = po.system

print("all subunits:",
      [a.details for a in s.asym_units])

print("first subunit sequence:",
      "".join(r.code for r in s.asym_units[0].entity.sequence))

print("all restraints on the system:", s.restraints)

import ihm.dumper
with open('initial.cif', 'w') as fh:
    ihm.dumper.write(fh, [s])

print("restraint datasets:", [r.dataset for r in s.restraints])

# Dataset for XL-MS restraint
d = s.restraints[0].dataset
print("XL-MS dataset at:", d.location.path)
print("Details:", d.location.details)

# Dataset for EM restraint
d = s.restraints[1].dataset
print("GMM file at", d.location.path)
print("is derived from EMDB entry", d.parents[0].location.access_code)

# Definitions of some common crosslinkers
import ihm.cross_linkers

# There should be exactly one XL restraint on the system
xl, = [r for r in s.restraints if isinstance(r, ihm.restraint.CrossLinkRestraint)]
xl.linker = ihm.cross_linkers.dss

# Get last step of last protocol (protocol is an 'orphan' because
# we haven't used it for a model yet)
last_step = s.orphan_protocols[-1].steps[-1]
print(last_step.num_models_end)

# Correct number of output models to account for multiple runs
last_step.num_models_end = 200000

import os

def get_cluster_members(analysis_dir, cluster_num):
    """Get filenames of RMF files in the given cluster"""
    num_models = 0
    for sample in ('A', 'B'):
        file_for_id = {}
        with open(os.path.join(analysis_dir,
                               'Identities_%s.txt' % sample)) as fh:
            for line in fh:
                filename, model_id = line.rstrip('\r\n').split()
                file_for_id[model_id] = filename
        with open(os.path.join(analysis_dir,
                               'cluster.%d.sample_%s.txt' % (cluster_num, sample))) as fh:
            for line in fh:
                model_id = line.rstrip('\r\n')
                num_models += 1
                yield os.path.join(analysis_dir, file_for_id[model_id])
                # In the interests of speed and space, let's just get
                # 10 models from each sample
                if num_models % 10 == 0:
                    break
all_cluster_0_models = list(get_cluster_members('../analysis', 0))

# Get last protocol in the file
protocol = s.orphan_protocols[-1]
# State that we filtered the 200000 frames down to one cluster:
import ihm.analysis
analysis = ihm.analysis.Analysis()
protocol.analyses.append(analysis)
analysis.steps.append(ihm.analysis.ClusterStep(
                      feature='RMSD', num_models_begin=200000,
                      num_models_end=len(all_cluster_0_models)))

mg = ihm.model.ModelGroup(name="Cluster 0")

# Add to last state
s.state_groups[-1][-1].append(mg)

# Make DCD of all models in cluster
dcd_filename = '../analysis/cluster.0/allmodels.dcd'
num_models = 0
with open(dcd_filename, 'wb') as dcd_fh:
    dcd = ihm.model.DCDWriter(dcd_fh)
    for model_file in all_cluster_0_models:
        rh = RMF.open_rmf_file_read_only(model_file)
        IMP.rmf.link_hierarchies(rh, [root_hier])
        IMP.rmf.load_frame(rh, RMF.FrameID(0))
        del rh
        m = po.add_model(mg)
        dcd.add_model(m)
        # We only want the model in DCD, not mmCIF
        del mg[-1]

dcd_location = ihm.location.OutputFileLocation(
    path=dcd_filename,
    details="Coordinates of all structures in the largest cluster")

e = ihm.model.Ensemble(model_group=mg,
                       num_models=len(all_cluster_0_models),
                       post_process=analysis.steps[-1],
                       name="Cluster 0", file=dcd_location)
s.ensembles.append(e)

# Add the cluster center model from RMF
rh = RMF.open_rmf_file_read_only('../analysis/cluster.0/cluster_center_model.rmf3')
IMP.rmf.link_hierarchies(rh, [root_hier])
IMP.rmf.load_frame(rh, RMF.FrameID(0))
del rh
m = po.add_model(mg)

em, = [r for r in s.restraints if isinstance(r, ihm.restraint.EM3DRestraint)]
em.fits[m] = ihm.restraint.EM3DRestraintFit(cross_correlation_coefficient=None)

for comp in ('C31', 'C34', 'C82', 'C53', 'C37'):
    asym = po.asym_units['%s.0' % comp]
    loc = ihm.location.OutputFileLocation('../analysis/cluster.0/LPD_%s.mrc' % comp)
    den = ihm.model.LocalizationDensity(file=loc, asym_unit=asym)
    e.densities.append(den)

datarepo = ihm.location.Repository(doi="10.5281/zenodo.3526621",
            root="../data", top_directory='data',
            url="https://zenodo.org/record/3526621/files/data.zip")
clusterrepo = ihm.location.Repository(doi="10.5281/zenodo.3526621",
            root="../analysis/cluster.0", top_directory='cluster.0',
            url="https://zenodo.org/record/3526621/files/cluster.0.zip")
s.update_locations_in_repositories([datarepo, clusterrepo])

# Dataset for XL-MS restraint
d = s.restraints[0].dataset
print("XL-MS dataset at %s/%s inside %s"
      % (d.location.repo.top_directory,
         d.location.path, d.location.repo.url))

# First localization density
d = e.densities[0]
print("Localization density at %s/%s inside %s"
       % (d.file.repo.top_directory,
          d.file.path, d.file.repo.url))

import ihm.dumper
with open('rnapoliii.cif', 'w') as fh:
    ihm.dumper.write(fh, [s])

#with open('rnapoliii.bcif', 'wb') as fh:
#    ihm.dumper.write(fh, [s], format='BCIF')

import ihm.reader
with open('rnapoliii.cif') as fh:
    s, = ihm.reader.read(fh)
print(s.title, s.restraints, s.ensembles, s.state_groups)

import urllib.request
with urllib.request.urlopen('https://pdb-dev.wwpdb.org/static/cif/PDBDEV_00000014.cif') as fh:
    s, = ihm.reader.read(fh)
print(s.title, s.restraints, s.ensembles, s.state_groups)
