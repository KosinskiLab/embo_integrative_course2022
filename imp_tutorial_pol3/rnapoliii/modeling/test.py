from __future__ import print_function

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
import pylab

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
if len(sys.argv)>1:
    output_prefix=A
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
                                   database=xl1,# The XLDB defined above
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


# ## Stage 3 - Sampling <a name="Sampling_2"></a>
# 
# With the system representation built and data restraints entered, the system is now ready to sample configurations. A replica exchange run can be set up using the `ReplicaExchange0` macro:
# 
# See the [ReplicaExchange0 documentation](https://integrativemodeling.org/nightly/doc/ref/classIMP_1_1pmi_1_1macros_1_1ReplicaExchange0.html) for a full description of all of the input parameters.
# 
# The sampling is performed by executing the macro:
# 
# ```mc1.execute_macro()```
# 

quick=False
if not quick:

    # total number of saved frames
    num_frames = 2

    # This object defines all components to be sampled as well as the sampling protocol
    mc1=IMP.pmi.macros.ReplicaExchange0(m,
              root_hier=root_hier,                         # The root hierarchy
              monte_carlo_sample_objects=dof.get_movers()+xl1rest.get_movers(), # All moving particles and parameters
              rmf_output_objects=outputobjects,            # Objects to put into the rmf file
              monte_carlo_temperature=1.0,   
              replica_exchange_minimum_temperature=1.0,
              replica_exchange_maximum_temperature=2.5,              
              simulated_annealing=False,
              number_of_best_scoring_models=0,
              monte_carlo_steps=10,
              number_of_frames=num_frames,
              save_coordinates_mode="25th_score",
              global_output_directory=output_directory)

    # Start Sampling
    mc1.execute_macro()


nmodels=10
rmffile="./results/150_xl_cryoem_3.rmf"

import IMP.pmi.output

hh=IMP.pmi.output.StatHierarchyHandler(m,rmffile)

#Total number of frames
print("Frames",len(hh))

# Describe the content of the first frame of the rmf file
print(hh[0])

#list down all the features names
for k in hh[0].features.keys(): print(k)
    

# For instance we can compute the distance between two residues


p0=IMP.atom.Selection(hh,molecule="C31",residue_index=10).get_selected_particles()[0]
p1=IMP.atom.Selection(hh,molecule="C34",residue_index=10).get_selected_particles()[0]

d0=IMP.core.XYZ(p0)
d1=IMP.core.XYZ(p1)

#note that hh can be used as a list
#pylab.plot([IMP.core.get_distance(d0,d1) for h in hh]);

#pylab.figure()

# Or we can get the radius of gyration of the whole complex
ps=IMP.atom.Selection(hh).get_selected_particles()
#pylab.plot([IMP.atom.get_radius_of_gyration(ps) for h in hh])

# To reduce I/O, we can store the data structure internal to hh, 
# so that it is not read directly from the files
# and it is faster

data=hh.data

# Then we plot the scores
#pylab.plot([x.score for x in data])

#pylab.figure() 

# finally we plot the distance of two crosslinked residues
#pylab.plot([float(x.features["CrossLinkingMassSpectrometryRestraint_Score_|XL|29.APO.1|C31|91|C160|1458|0|CLASS_0|"]) for x in data]);


scores={}

for xl in xl1: 
    if not xl['IntraRigidBody']:
        scores.update({xl['XLUniqueSubID']:float(xl['IDScore'])})
    

x=[]
y=[]
for k in data[0].features.keys():
    if "Distance" in k:
        id=k.split("|")[2]
        if id in scores:
            x.append(scores[id])
            y.append(float(data[2].features[k]))

#pylab.scatter(x,y);


are=IMP.pmi.macros.AnalysisReplicaExchange(m,
                 [rmffile],
                 best_models=nmodels,
                 alignment=False)

print(are)


are.set_rmsd_selection(molecules=["C31","C34","C53","C37","C82"])

are.cluster(10.0)

# see the contant of the "are" object
print(are)

#print the cluster info
for cluster in are:
    print(cluster)

for member in are[0]:
    print(member)
    
are.save_coordinates(are[0])

# slow!
# are.plot_rmsd_matrix("rmsd_matrix.pdf")


keys=[k for k in are[0][0].features.keys() if "CrossLinkingMassSpectrometryRestraint" in k]
#pylab.figure()
for ind in ["0","1","2"]:
    distances={}
    psi=[]
    sigma=[]
    for member in are[0]:
        for k in keys:
            if "CLASS_"+str(ind) in k and "Distance" in k:
                if k in distances:
                    distances[k].append(float(member.features[k]))
                else:
                    distances[k]=[float(member.features[k])]
            if "CLASS_"+str(ind) in k and "Psi" in k and "MonteCarlo" not in k:
                psi.append(float(member.features[k]))
            if "SIGMA" in k and "Sigma" in k and "MonteCarlo" not in k:
                sigma.append(float(member.features[k]))

    x=[distances[k] for k in distances]
    #pylab.boxplot(x);
    #pylab.figure()
    #pylab.plot(range(len(psi)),psi);
    #pylab.figure()
#pylab.plot(range(len(sigma)),sigma);




### build the seed 

are.write_seed("seed.rmf3", 2)

are.compute_cluster_center(cluster=are[0])
print(are.precision(cluster=are[0]))
print(are.bipartite_precision(cluster1=are[0],cluster2=are[1]))

rmsf1=are.rmsf(cluster=are[0],molecule='C31');
#pylab.plot(list(rmsf1.keys()),list(rmsf1.values()),marker=".",linewidth=0)
#pylab.figure()

rmsf2=are.rmsf(cluster=are[0],molecule='C34');
#pylab.plot(list(rmsf2.keys()),list(rmsf2.values()),marker=".",linewidth=0)

for mol in ['ABC23','ABC10beta','ABC14_5','ABC27','C25','AC40','C160','ABC10alpha','C128','AC19','C11','C17','C31','C34','C53','C37','C82']: 
    are.rmsf(cluster=are[0],molecule=mol);
ch1=IMP.pmi.tools.ColorHierarchy(are.stath1)
ch1.color_by_uncertainty()
are.save_coordinates(are[0])

density_names={'core': ['ABC23','ABC10beta','ABC14_5','ABC27','C25','AC40','C160','ABC10alpha','C128','AC19','C11','C17'],
               'C53': ['C53'], 
               'C37': ['C37'], 
               'C34': ['C34'], 
               'C82': ['C82'], 
               'C31': ['C31']}

# you can iterate on the clusters
for n,a in enumerate(are):
    are.save_densities(cluster=a,density_custom_ranges=density_names,prefix="Cluster-"+str(n))

# it is slow
# are.contact_map(cluster=are[0]);


