#----------------------------------------------------------------------------------------
# This script is tetsed for Chimera 1.11.2
# Called as:
#
#    chimera --nogui --script "CCcalculate <pdb file> <map> <resolution> <n searches>"
#
# A synthetic map of the PDB structure with the resolution of the experimental map
# will be generated and fitted on the experimental map with n searches.
# After fitting, the correlation of the best fitted structure is measured
#----------------------------------------------------------------------------------------

import argparse
import chimera
import VolumeViewer
from MoleculeMap import molmap
from FitMap.fitcmd import fitmap
from Measure import measure

# Parsing information
parser = argparse.ArgumentParser(description="Fits model in EM density map using n gloval fit searches")
parser.add_argument("model", help="Path to model PDB")
parser.add_argument("em_map", help="Path to experimental EM density map")
parser.add_argument("resolution", help="Resolution of EM map")
parser.add_argument("searches", help="Number of global fit searches")
args = parser.parse_args()

model = chimera.openModels.open(args.model)[0]
em_map = VolumeViewer.open_volume_file(args.em_map)[0]
res = float(args.resolution)
search = int(args.searches)

# Synthetic map generation and global fitting

synth_map = molmap.molecule_map(model.atoms, res)
select = chimera.selection.ItemizedSelection([synth_map])
fit_list = fitmap(select, em_map, search = search, listFits = False)

measure.correlation("correlation", [synth_map], [em_map])
