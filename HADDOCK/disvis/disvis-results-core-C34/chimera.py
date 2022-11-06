from chimera import runCommand
from chimera import openModels
import Midas
import os
print(os.getcwd())

sides= {    'front'  :[ 1, 0, 0,   0], \
            'right'  :[ 0, 1, 0,  90], \
            'back'   :[ 0, 1, 0, 180], \
            'left'   :[ 0, 1, 0, 270], \
            'bottom' :[ 1, 0, 0,  90], \
            'top'    :[ 1, 0, 0, 270] }

#set bg_color white
runCommand('preset apply pub 3')
openModels.open('fixed_chain.pdb')
vol=openModels.open('accessible_interaction_space.mrc')[0]

runCommand('color dark slate gray ,s ; surftransparency 50')
runCommand('sop smooth # factor 0.5 iterations 3 inPlace true')
runCommand('color firebrick,r helix; color goldenrod,r strand; color gray,r coil')

Midas.center('#')
Midas.window('#')

for i in range(1,int(vol.matrix_stats.maximum)+1):
    vol.set_parameters(surface_levels=[i])
    vol.show()
    for side in sides.keys():
        Midas.turn(sides[side][:3], angle=sides[side][3])

        Midas.center('#')
        Midas.window('#')
        Midas.copy(file="{}_{}.png".format(i, side), width=400, height=400)

        Midas.turn(sides[side][:3], angle=-sides[side][3])
quit()

#runCommand('center; window; scale 0.95')

#runCommand('copy file img1.png width 800 height 800 ')
#runCommand('turn 1,0,0 90 ')
#runCommand('center; window; scale 0.95')
#runCommand('copy file img2.png width 400 height 400')
#runCommand('turn 0,1,0 90 ')
#runCommand('center; window; scale 0.95')
#runCommand('copy file img3.png width 800 height 800 ')
#runCommand('turn 0,0,1 90 ')
#runCommand('center; window; scale 0.95')
#runCommand('copy file img4.png width 800 height 800 ')
#runCommand('stop')
