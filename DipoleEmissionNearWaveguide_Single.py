import meep as mp
from meep import mpb
import numpy as np
from matplotlib import pyplot as plt

omega = 1/0.78
# index of material
index_SiN=2.02
index_SiO2=1.46
index_PVA=1.5
index_An=1.8
# The size of the guide
width=0.65
height=0.25

height_An=0.2
height_PVA=0.03

# The size of simulation domain
width_cal=2.3
length_cal=2
height_cal=2.5

# Thinkness of PML
t_pml=0.4

# The size of box to obatin the dipole power
size_x=0.3
size_y=0.3
size_z=0.3

# The position of dipole
pos_x=0
pos_y=0
pos_z=height+0.005

pt= mp.Vector3( pos_x, pos_y, pos_z)

resolution=200
numx=int(width_cal*resolution)
numy=int(height_cal*resolution)
x_list=np.linspace(-width_cal/2,width_cal/2,numx)
y_list=np.linspace(-height_cal/2,height_cal/2,numy)
x_grid,y_grid=np.meshgrid(x_list,y_list)
# %%
cell_size = mp.Vector3(length_cal,width_cal,height_cal)

pml_layers = [mp.Absorber(thickness=t_pml,direction=mp.Y),
                   mp.PML(thickness=t_pml,direction=mp.Z),
                   mp.PML(thickness=t_pml,direction=mp.X)]

# Y dipole 
source = [mp.Source(mp.ContinuousSource(frequency= omega), 
            component = mp.Ey, center=pt )]

geometry = [mp.Block(material=mp.Medium(epsilon=index_SiN**2), 
                     size=mp.Vector3(mp.inf,width, height),
                     center=mp.Vector3(0, 0, height/2)),

            mp.Block(material=mp.Medium(epsilon=index_SiO2**2),
                     size=mp.Vector3(mp.inf,mp.inf, height_cal/2),
                     center=mp.Vector3(0, 0, -height_cal/4)),

            mp.Block(material=mp.Medium(epsilon=index_An**2),
                     size=mp.Vector3(mp.inf,mp.inf, height_An),
                     center=mp.Vector3(0, 0, height+height_An/2)),
            mp.Block(material=mp.Medium(epsilon=index_PVA**2),
                     size=mp.Vector3(mp.inf,mp.inf, height_PVA),
                     center=mp.Vector3(0, 0, height+height_An+height_PVA/2))
                     ]

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=source,
                    eps_averaging= False,
                    Courant = 0.5)

box_1 = sim.add_mode_monitor(omega, 0, 1,
mp.FluxRegion(center=mp.Vector3(0.2, 0, 0),
                         size=mp.Vector3(0,width_cal,height_cal)))  
box_2 = sim.add_mode_monitor(omega, 0, 1,
mp.FluxRegion(center=mp.Vector3(-0.2, 0, 0),
                         size=mp.Vector3(0,width_cal,height_cal))) 


box_x1 = sim.add_flux(omega, 0, 1,
mp.FluxRegion(center=mp.Vector3(pos_x-size_x/2,pos_y,pos_z),size=mp.Vector3(0,size_y,size_z)))
box_x2 = sim.add_flux(omega, 0, 1,
mp.FluxRegion(center=mp.Vector3(pos_x+size_x/2,pos_y,pos_z),size=mp.Vector3(0,size_y,size_z)))

box_y1 = sim.add_flux(omega, 0, 1,
mp.FluxRegion(center=mp.Vector3(pos_x,pos_y-size_y/2,pos_z),size=mp.Vector3(size_x,0,size_z)))
box_y2 = sim.add_flux(omega, 0, 1,
mp.FluxRegion(center=mp.Vector3(pos_x,pos_y+size_y/2,pos_z),size=mp.Vector3(size_x,0,size_z)))

box_z1 = sim.add_flux(omega, 0, 1,
mp.FluxRegion(center=mp.Vector3(pos_x,pos_y,pos_z-size_z/2),size=mp.Vector3(size_x,size_y,0)))
box_z2 = sim.add_flux(omega, 0, 1,
mp.FluxRegion(center=mp.Vector3(pos_x,pos_y,pos_z+size_z/2),size=mp.Vector3(size_x,size_y,0)))

# %%
# sim.run(until_after_sources=
# mp.stop_when_fields_decayed(20, mp.Ey, pt, 1e-3))
# sim.run(until_after_sources=20)
sim.run(until=20)
# %%
x1 = mp.get_fluxes(box_x1)
x2 = mp.get_fluxes(box_x2)
y1 = mp.get_fluxes(box_y1)
y2 = mp.get_fluxes(box_y2)
z1 = mp.get_fluxes(box_z1)
z2 = mp.get_fluxes(box_z2)

mode1 = sim.get_eigenmode_coefficients(box_1,
                      [1], eig_parity=mp.NO_PARITY)
mode2 = sim.get_eigenmode_coefficients(box_2,
                      [1], eig_parity=mp.NO_PARITY)

# %%
ptotal = -x1[0]+x2[0]-y1[0]+y2[0]-z1[0]+z2[0]
pmode_1 = abs(mode1.alpha[0,0,0])**2
pmode_2 = abs(mode2.alpha[0,0,1])**2
effic_1=pmode_1/ptotal
effic_2=pmode_2/ptotal
neff1 = mode1.kdom[0].x/omega
neff2 = mode2.kdom[0].x/omega
print(neff1)
print(neff2)
print("efficency to position direction:, {} ".format(effic_1))
print("efficency to negative direction:, {} ".format(effic_2))
# %%
