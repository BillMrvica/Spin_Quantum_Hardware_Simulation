# due to the way KLayout handles Python imports, to import module patching.py
# we firstly have to add the absolute current working directory to sys.path
import os
import sys
script_path = os.path.abspath(__file__)
current_dir = os.path.dirname(script_path)
if current_dir not in sys.path:
    sys.path.append(current_dir)

import pya

# -----------------------------------------------------------------------------
from importlib import reload
import PatchingDialog
reload(PatchingDialog)  # reload the module to ensure it's up-to-date

from PatchingDialog import PatchingDialog

dialog = PatchingDialog() # Create an instance of dialog
dialog.exec_() # Show the dialog and awaits until the user closes it


# # --------------------------- Debug mode ---------------------------
# from importlib import reload
# import patching
# reload(patching)  # reload the module to ensure it's up-to-date

# from patching import create_grid_from_shapes, create_patch

# # standard code to get the current layout and cell: 
# # application instance -> main window -> current view -> active cellview -> layout and active cell
# active_cellview = pya.Application.instance().main_window().current_view().active_cellview()
# layout = active_cellview.layout()
# top_cell = active_cellview.cell # by default, you should be working on the top cell

# # user input
# electrode_layer = 111
# box_layer = 9
# patch_layer = 113
# patch_size = 5.0 # um

# # create grid
# grid = create_grid_from_shapes(layout=layout,
#                                 cell=top_cell,
#                                 shape_layer=box_layer,  # BOX_LAYER/0 is the layer where boxes are drawn
#                                 grid_line_width=patch_size,)  # width of the grid lines in micrometers

# # find intersections and create patches
# create_patch(layout=layout,
#               cell=top_cell,
#               grid_region=grid,
#               electrode_layer=electrode_layer,  # ELECTRODE_LAYER/0 is the layer where electrodes are drawn
#               patch_layer=patch_layer,  # PATCH_LAYER/0 is the layer where patches will be created
#               patch_size=patch_size)  # size of each patch in micrometers
