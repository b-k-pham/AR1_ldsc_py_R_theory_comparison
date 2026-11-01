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

# python/ (this file's directory) contains src/; data/ is a sibling of python/ at the repo root
PROJECT_ROOT = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

DATA_DIR = PROJECT_ROOT.parent / "data"

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


### MAKE REF PANEL
n = 2000
X_properties = {'n': n,'m':n*2,'pm_causal':None,'fixed_m_causal': False,'sigma_s':0, 'rho1':0.3,'rho2':None,'Fst':0}
ref_X_properties = {'ref_n': n,'ref_pm_causal':None,'ref_fixed_m_causal': False,'ref_sigma_s':0.3, 'ref_rho1':0.3,'ref_rho2':None,'ref_Fst':0.01}
ld_mat_properties = {'realistic': False,'model_Fst_in_realistic':True,'prefix': '1kg_p1_eur_ben/1kg_p1_eur_chr1','make_ref_ldscores':True}
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

#rhos = [0.995]
counter = -1
multithreading = True
n_jobs = 1
sigma_s = 0.3
Fst = np.float64(sys.argv[1])
Fst_label = str(Fst).replace(".", "") 

ref_data_gen = data_generation(n = ref_n,m= m,Fst = Fst,rho1 = ref_rho1,rho2 = ref_rho2,sigma_s = sigma_s,h2_pop = h2_pop,pm_causal = pm_causal,fixed_m_causal = ref_fixed_m_causal, old = old,realistic = realistic,prefix = prefix,use_gwash_m = False,do_ldscore_bias_correct = True)
ref_ldscores,ref_mu2_hat,ref_mu3_hat,X_ref,f_ref = gen_ref_ldscores_panel(ref_data_gen,seed = seed_num,nPCs = nPCs, regress_X_on_PC = regress_X_on_PC, regress_y_on_PC = regress_y_on_PC)
n_tilde = ref_data_gen.n
# with open('ref_panel_Fst0005_sigma_s_06_chr1.npy', 'wb') as f:
#     np.save(f,{'ref_ldscores':ref_ldscores,'ref_mu2_hat':ref_mu2_hat,'ref_mu3_hat':ref_mu3_hat,'X_ref':X_ref})

print(np.mean(ref_ldscores))

with open(DATA_DIR / "ref_panel_npy" / 'n{n}_m{m}_ref_panel_Fst{Fst_label}_sigma_s_03_rho03_AR1_sstats_only_seed123.npy'.format(n = n,m = m,Fst_label = Fst_label), 'wb') as f:
    np.save(f,{'ref_data_gen': ref_data_gen,'ref_ldscores':ref_ldscores,'ref_mu2_hat':ref_mu2_hat,'ref_mu3_hat':ref_mu3_hat,'f_ref':f_ref})