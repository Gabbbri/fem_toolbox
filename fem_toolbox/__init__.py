# EXPOSING SUBMODULES 


from . import geometry
from . import femsolver
from . import BC_loads
from . import elements
from . import mesh
from . import postprocessing
from . import optimizer


from pyfiglet import Figlet

fig = Figlet(font='doom')
print(fig.renderText('fem_toolbox'))
print("github.com/Gabbbri/fem_toolbox")

