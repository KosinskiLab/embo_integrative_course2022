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

m=IMP.Model()

nmodels=150
rmffile="150_xl_cryoem_3.rmf"

import IMP.pmi.output


are=IMP.pmi.macros.AnalysisReplicaExchange(m,
                 [rmffile],
                 best_models=nmodels,
                 alignment=True)

print(are)


are.set_rmsd_selection(molecules=["C31","C34","C53","C37","C82"])
are.set_alignment_selection(molecules=['ABC23','ABC10beta','ABC14_5','ABC27','C25','AC40','C160','ABC10alpha','C128','AC19','C11','C17']+["C31","C34","C53","C37","C82"])

#are.cluster(20.0)

# see the contant of the "are" object
#print(are)

#print the cluster info
#maxcluster=None
#maxsize=0
#for cluster in are:
#    if len(cluster)>maxsize:
#        maxcluster=cluster


are[0].center_index=0


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
are.save_densities(cluster=are[0],density_custom_ranges=density_names,prefix="Ensemble")



