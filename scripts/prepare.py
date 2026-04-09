import os
from calvados.cfg import Config, Job, Components
import subprocess
import numpy as np
from argparse import ArgumentParser
from Bio import SeqIO

parser = ArgumentParser()
parser.add_argument('--name',nargs='?',required=True,type=str)
parser.add_argument('--ionic',nargs='?',required=True,type=int)
args = parser.parse_args()

cwd = os.getcwd()
N_save = int(5e4)
N_frames = 9000

sysname = f'{args.name:s}_{args.ionic:d}'
residues_file = f'{cwd}/input/residues_CALVADOS2.csv'

config = Config(
  # GENERAL
  sysname = sysname, # name of simulation system
  box = [20, 20, 300.], # nm
  temp = 298.15,
  ionic = args.ionic/1000, # molar
  pH = 7.5,
  topol = 'slab',
  slab_width = 30,
  friction = 0.01,

  # RUNTIME SETTINGS
  wfreq = N_save, # dcd writing frequency, 1 = 10fs
  steps = N_frames*N_save, # number of simulation steps
  runtime = 0, # overwrites 'steps' keyword if > 0
  platform = 'CUDA' if N_frames > 0 else 'CPU',
  restart = 'checkpoint',
  frestart = 'restart.chk',
  verbose = True,
  slab_eq = True,
  steps_eq = 100*N_save,
)

# PATH
path = f'{cwd}/{sysname}'
output_path = f'{cwd}/data'
subprocess.run(f'mkdir -p {path}',shell=True)
subprocess.run(f'mkdir -p {output_path}',shell=True)

analyses = f"""
from calvados.analysis import SlabAnalysis, calc_conf_profiles, calc_com_traj, calc_contact_map

slab = SlabAnalysis(name="{sysname:s}", input_path="{path:s}",
                    output_path="{output_path:s}", ref_name="{sysname:s}", verbose=True)

slab.center(start=400, center_target='all')
slab.calc_profiles()
slab.calc_concentrations()
calc_com_traj(path="{path:s}",sysname="{sysname:s}",output_path="{output_path:s}",
        residues_file="{residues_file:s}",start=400)
calc_contact_map(path="{path:s}",sysname="{sysname:s}",output_path="{output_path:s}",is_slab=True)
"""

config.write(path,name='config.yaml',analyses=analyses)

components = Components(
  # Defaults
  molecule_type = 'protein',
  nmol = 1, # number of molecules
  restraint = False, # apply restraints
  charge_termini = 'both', # charge N or C or both

  # INPUT
  ffasta = f'{cwd}/input/fastalib.fasta', # input fasta file
  fresidues = residues_file, # residue definitions
)

components.add(name=args.name, nmol=200)

components.write(path,name='components.yaml')

