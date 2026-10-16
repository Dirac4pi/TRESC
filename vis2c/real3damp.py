#!/home/lky/miniconda3/envs/vis2c/bin/python
'''
Visualize a real-space scalar field, with optional cube-based color mapping.
Note: This is capable of plotting a non-axis-aligned cube.
author:Dirac4pi
env:vis2c
'''

import argparse
from os import path

import numpy as np
import viskit as vk
import plot3d as p3
from mayavi import mlab
from traits.api import HasTraits, Float, observe
from traitsui.api import View, Item

#-------------------------------------------------------------------------------
def real3damp(ampcube: str, isovalue: float = 0.05) -> None:
  """Plot one cube field without color mapping.

  Parameters:
    ampcube: Gaussian .cub/.cube file containing the real-space scalar field.
    isovalue: Positive level for the single isosurface. The plotting
      function may adjust it according to its current threshold rules.

  Returns:
    None. Opens a TraitsUI controller for changing the isovalue.
  """
  if not np.isfinite(isovalue):
    raise ValueError("isovalue must be finite.")
  if path.splitext(ampcube)[1].lower() not in (".cub", ".cube"):
    raise ValueError("ampcube must be a .cub or .cube file.")
  ampbase = path.splitext(path.basename(ampcube))[0]
  print(f"Loading real-space field from {ampcube}...")
  atoms, ampdat = vk.load_cube(ampcube)
  x, y, z = ampdat["x"], ampdat["y"], ampdat["z"]
  values = ampdat["value"]
  if not np.all(np.isfinite(values)):
    raise ValueError("Cube contains non-finite field values.")
  # This title triggers atom plotting inside real_func_ampplot_cub().
  title = f"amplitude: {ampbase}"
  def plot(level):
    """Draw one isosurface using the same function on every update."""
    used_level = p3.real_func_ampplot_cub(
      title, atoms, values, x, y, z, level
    )
    return used_level, mlab.figure(title)
  print("Plotting...")
  used_level, figure = plot(isovalue)
  mlab.draw(figure=figure)

  class IsoValueController(HasTraits):
    isov = Float(used_level, desc="isov", auto_set=False, enter_set=True)
    @observe("isov")
    def update_isov(self, event):
      new_level = event.new
      if not np.isfinite(new_level):
        print("isovalue must be finite.")
        return
      try:
        _, updated_figure = plot(new_level)
      except ValueError as exc:
        print(f"Cannot update isovalue: {exc}")
        return
      mlab.draw(figure=updated_figure)
    view = View(
      Item("isov", label="isovalue", show_label=True),
      width=300,
      height=200,
      resizable=True
    )

  controller = IsoValueController()
  controller.configure_traits()
  return None

#===============================================================================
def main(argv=None) -> None:
  """Parse CLI arguments and launch a single-cube display."""
  parser = argparse.ArgumentParser(
    description='Plot a cube isosurface of a real-space function'
  )
  parser.add_argument('ampcube', help='Amplitude cube (.cub/.cube).')
  parser.add_argument('--isovalue', type=float, default=0.05,
                      help='Signed contour level (default: 0.05).')
  args = parser.parse_args(argv)
  try:
    real3damp(args.ampcube, args.isovalue)
  except (ValueError, FileNotFoundError, NotImplementedError) as exc:
    parser.exit(2, f'{parser.prog}: {exc}\n')
  return None


if __name__ == '__main__':
  main()
