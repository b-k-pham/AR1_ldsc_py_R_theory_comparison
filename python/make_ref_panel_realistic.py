import pandas as pd
import numpy as np
import scipy
import scipy.sparse as sp
import scipy.io as sio
import scipy.stats as stats
from tqdm.notebook import tqdm


import os
os.environ["R_HOME"] = f"{os.environ['CONDA_PREFIX']}\\Lib\\R"


from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

from plotnine import *

import matplotlib.pyplot as plt

import pickle

from joblib import Parallel, delayed

import sys
from pathlib import Path

# Block dask's optional `distributed` import: in this env, importing it pulls in tornado,
# whose SSL cert loading raises (not ImportError), which dask's except-ImportError guard
# doesn't catch. Pre-poisoning sys.modules makes the import raise ImportError instead.
sys.modules['distributed'] = None

# python/ (this file's directory) contains src/; data/ and 1kg_p1_eur_ben/ are siblings of
# python/ at the repo root.
PROJECT_ROOT = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

REPO_ROOT = PROJECT_ROOT.parent
DATA_DIR = REPO_ROOT / "data"
PLINK_PREFIX = str(REPO_ROOT / "1kg_p1_eur_ben" / "1kg_p1_eur_chr1")

import dask
dask.config.set(scheduler='synchronous')  # avoid importing distributed/tornado (broken SSL cert load in this env)

from src.methods.GWASH_funcs import *
from src.methods.GWASH_sim_funcs import *
from src.methods.ldsc_barebones import *

from src.simtools.sim_utils import *
from src.simtools.data_generation import data_generation as data_generation
from src.simtools.data_preprocessing import data_preprocessing as data_preprocessing
from src.simtools.do_analysis import do_analysis as do_analysis
from src.simtools.run_simulations import run_simulations as run_simulations
from src.simtools.data_generation import gen_ref_ldscores_panel as gen_ref_ldscores_panel
from src.simtools.visualization import visualize

from natsort import natsorted

import pickle

seed_num = 123
np.random.seed(seed_num)


### MAKE REF PANEL, LD STRUCTURE TAKEN FROM REAL 1000G chr1 GENOTYPES (realistic=True)
# instead of synthetic AR1. Everything else mirrors make_ref_panel.py.
# Usage: python make_ref_panel_realistic.py <Fst> [n]   (m is always n*2; default n=500)
n = int(sys.argv[2]) if len(sys.argv) > 2 else 500
X_properties = {'n': n,'m':n*2,'pm_causal':None,'fixed_m_causal': False,'sigma_s':0, 'rho1':0.3,'rho2':None,'Fst':0}
ref_X_properties = {'ref_n': n,'ref_pm_causal':None,'ref_fixed_m_causal': False,'ref_sigma_s':0.3, 'ref_rho1':0.3,'ref_rho2':None,'ref_Fst':0.01}
ld_mat_properties = {'realistic': True,'model_Fst_in_realistic':True,'prefix': PLINK_PREFIX,'make_ref_ldscores':True}
method_properties = {'reml_tol': 1e-8,'reml_max_iters': 100}
stats_properties = {'scaleX': True,'scaley': True, 'nPCs':5,'regress_PC_out':False,'regress_X_on_PC': False,'regress_y_on_PC': False}
simul_properties = {'num_sims':np.nan,'h2_pop': 0.3}
debug_properties = {'calc_mu_hat_2_fast': True,'old':False,'track_progress':True,'run_gcta':False}

list_of_dicts = [X_properties,ref_X_properties,ld_mat_properties,method_properties,stats_properties,simul_properties,debug_properties]
params  = combine_all_dicts(list_of_dicts)

to_run = dict()
to_run['demonstration'] = params

sim_key = 'demonstration'
locals().update(to_run[sim_key])
res_dict = dict()
res_dict_raw = dict()

counter = -1
multithreading = True
n_jobs = 1
sigma_s = 0.3
Fst = np.float64(sys.argv[1])
Fst_label = str(Fst).replace(".", "")

ref_data_gen = data_generation(n = ref_n,m= m,Fst = Fst,rho1 = ref_rho1,rho2 = ref_rho2,sigma_s = sigma_s,h2_pop = h2_pop,pm_causal = pm_causal,fixed_m_causal = ref_fixed_m_causal, old = old,realistic = realistic,prefix = prefix,use_gwash_m = False,do_ldscore_bias_correct = True)
ref_ldscores,ref_mu2_hat,ref_mu3_hat,X_ref,f_ref = gen_ref_ldscores_panel(ref_data_gen,seed = seed_num,nPCs = nPCs, regress_X_on_PC = regress_X_on_PC, regress_y_on_PC = regress_y_on_PC)
n_tilde = ref_data_gen.n

print(np.mean(ref_ldscores))

with open(DATA_DIR / "ref_panel_npy_realistic" / 'n{n}_m{m}_ref_panel_Fst{Fst_label}_sigma_s_03_realisticLD_chr1_seed123.npy'.format(n = n,m = m,Fst_label = Fst_label), 'wb') as f:
    np.save(f,{'ref_data_gen': ref_data_gen,'ref_ldscores':ref_ldscores,'ref_mu2_hat':ref_mu2_hat,'ref_mu3_hat':ref_mu3_hat,'f_ref':f_ref})
