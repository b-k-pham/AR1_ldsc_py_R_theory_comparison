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

n = 2000
seed_num = 123
np.random.seed(seed_num)

X_properties = {'n': n,'m':n*2,'pm_causal':None,'fixed_m_causal': False,'sigma_s':0, 'rho1':0.3,'rho2':None,'Fst':0}
ref_X_properties = {'ref_n': n,'ref_pm_causal':None,'ref_fixed_m_causal': False,'ref_sigma_s':0, 'ref_rho1':0.3,'ref_rho2':None,'ref_Fst':0}
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

res_dfs = []

Fst = np.float64(sys.argv[1])
Fst_label = str(Fst).replace(".", "")
for i in tqdm(range(100)):
    if make_ref_ldscores:
        ref_dat = np.load(DATA_DIR / "ref_panel_npy" / 'n{n}_m{m}_ref_panel_Fst{Fst_label}_sigma_s_03_rho03_AR1_sstats_only_seed123.npy'.format(n = n,m = m,Fst_label = Fst_label),allow_pickle = True).item()
        ref_data_gen,ref_ldscores,ref_mu2_hat,ref_mu3_hat,f_ref = ref_dat['ref_data_gen'],ref_dat['ref_ldscores'],ref_dat['ref_mu2_hat'],ref_dat['ref_mu3_hat'],ref_dat['f_ref']
        n_tilde = ref_data_gen.n
    else:
        ref_data_gen = None
        ref_ldscores = None
        n_tilde = None
        ref_mu2_hat = None
        ref_mu3_hat = None
        X_ref = None
        f_ref = None

    # INDEPENDENT-F VARIANT: unlike do_AR1_sim_Sam_ref.py, we deliberately do NOT pass the
    # reference panel's f_ref into data_generation - ref_f stays None, so gen_X draws its own
    # fresh, independent f for every one of the 100 replicates (matching the pre-fix behavior).
    # This reproduces the "reference panel vs. in-sample f mismatch" scenario for the report's
    # second section, without touching the current with_ref_panel/ (fixed-f) data.
    my_data_gen = data_generation(n = n,m= m,Fst = Fst,rho1 = rho1,rho2 = rho2,sigma_s = sigma_s,h2_pop = h2_pop,pm_causal = pm_causal,fixed_m_causal = fixed_m_causal, old = old,model_Fst_in_realistic = model_Fst_in_realistic,realistic = realistic,prefix = prefix,ref_ldscores = ref_ldscores,ref_mu2_hat = ref_mu2_hat,ref_mu3_hat = ref_mu3_hat,ref_f = None)
    np.random.seed(seed_num + i)

    X,b,f,_ = my_data_gen.gen_X(realistic = realistic,debug = True)
    Cf = np.mean(f**2/(1+f**2))
    y = my_data_gen.gen_y(X,b)


    my_data_preprocessing = data_preprocessing(my_data_gen,nPCs = nPCs,regress_X_on_PC = regress_X_on_PC,regress_y_on_PC = regress_y_on_PC)
    if regress_PC_out:
        regress_res = my_data_preprocessing.regress_out_PCs(X,y)
        X = regress_res['X_res']
        y = regress_res['y_res']

    h2_samp_fve =  np.var(np.matmul(X,b),ddof=1)/np.var(y,ddof = 1)
    if scaleX:
        X_tilde = my_data_preprocessing.manual_standardize_and_scale_inR(X)
    else:
        X_tilde = X
    if scaley:
        y_tilde = my_data_preprocessing.manual_standardize_and_scale_inR(y)
    else:
        y_tilde = y
    # if my_data_gen.ref_ldscores is None, then generate ldscores. Otherwise, use my_data_gen.ref_ldscores
    u2,s2,ldscores = my_data_gen.compute_u2_ldscores(X_tilde,y_tilde)
    ldsc_reg_weights = ldscores
    sims = do_analysis(my_data_gen,X_tilde,u2,s2,ldscores,ldsc_reg_weights,calc_mu_hat_2_fast = True)
    sims_fixed = do_analysis(my_data_gen,X_tilde,u2,s2,ldscores,ldsc_reg_weights,intercept = 1, calc_mu_hat_2_fast = True)

    # VANILLA LDSC
    h2_ldsc_reg,se_ldsc_reg,icpt_ldsc_reg,weights = sims.do_ldsc()
    ldsc_reg_res = pd.DataFrame({"h2_ldsc_reg": h2_ldsc_reg, "se_ldsc_reg": se_ldsc_reg, "icpt_ldsc_reg": icpt_ldsc_reg},index = [i])
    # LDSC FIXED
    h2_ldsc_fixed,se_ldsc_fixed,icpt_ldsc_fixed,weights_fixed = sims_fixed.do_ldsc()
    ldsc_fixed_res = pd.DataFrame({"h2_ldsc_fixed": h2_ldsc_fixed, "se_ldsc_fixed": se_ldsc_fixed, "icpt_ldsc_fixed": icpt_ldsc_fixed},index = [i])

    # NO WEIGHTING LDSC FIXED
    ldsc_weights_ones = np.ones(m).reshape(-1,1)
    sims_no_weighting = do_analysis(my_data_gen,X_tilde,u2,s2,ldscores,ldsc_weights_ones,intercept = 1,num_weight_iteration = 0,calc_mu_hat_2_fast = True,no_weights = True)
    h2_ldsc_fixed_nw,se_ldsc_fixed_nw,icpt_ldsc_fixed_nw,weights_fixed_nw = sims_no_weighting.do_ldsc()
    ldsc_fixed_nw_res = pd.DataFrame({"h2_ldsc_fixed_nw": h2_ldsc_fixed_nw, "se_ldsc_fixed_nw": se_ldsc_fixed_nw, "icpt_ldsc_fixed_nw": icpt_ldsc_fixed_nw},index = [i])

    # NO WEIGHTING IN LDSC (NUMBER OF WEIGHTING ITERATIONS = 0 (default 2) and supplied weights are all 1 for each SNP
    ldsc_weights_ones = np.ones(m).reshape(-1,1)
    sims_no_weighting = do_analysis(my_data_gen,X_tilde,u2,s2,ldscores,ldsc_weights_ones,num_weight_iteration = 0,calc_mu_hat_2_fast = True,no_weights = True)
    h2_ldsc_nw,se_ldsc_nw,icpt_ldsc_nw,weights_ldsc_nw = sims_no_weighting.do_ldsc()
    ldsc_nw_res = pd.DataFrame({"h2_ldsc_nw": h2_ldsc_nw, "se_ldsc_nw": se_ldsc_nw, "icpt_ldsc_nw": icpt_ldsc_nw},index = [i])

    # NO TWOSTEP IN LDSC (TWOSTEP = np.inf (default = 30) SO EVERY SNP IS COUNTED) and stop at step 1.
    sims_no_ts = do_analysis(my_data_gen,X_tilde,u2,s2,ldscores,ldsc_reg_weights,calc_mu_hat_2_fast = True,twostep = np.inf,return_step1_only = True)
    h2_ldsc_nts,se_ldsc_nts,icpt_ldsc_nts,weights_ldsc_nts = sims_no_ts.do_ldsc()
    ldsc_nts_res = pd.DataFrame({"h2_ldsc_nts": h2_ldsc_nts, "se_ldsc_nts": se_ldsc_nts, "icpt_ldsc_nts": icpt_ldsc_nts},index = [i])

    # NO WEIGHTING AND NO TWOSTEP IN LDSC
    sims_no_weighting_and_no_ts = do_analysis(my_data_gen,X_tilde,u2,s2,ldscores,ldsc_weights_ones,num_weight_iteration = 0,calc_mu_hat_2_fast = True,no_weights = True,twostep = np.inf,return_step1_only = True)
    h2_ldsc_nwts,se_ldsc_nwts,icpt_ldsc_nwts,weights_ldsc_nwts = sims_no_weighting_and_no_ts.do_ldsc()
    ldsc_nwts_res = pd.DataFrame({"h2_ldsc_nwts": h2_ldsc_nwts, "se_ldsc_nwts": se_ldsc_nwts, "icpt_ldsc_nwts": icpt_ldsc_nwts},index = [i])

    if make_ref_ldscores is False:
        gwash_res_sample_theoretical = sims.do_GWASH_from_X_tilde() # Do same calculation from X_tilde() just to get theoretical variance.
        gwash_mu2_hat = gwash_res_sample_theoretical['mu2_hat']
        gwash_mu3_hat = gwash_res_sample_theoretical['mu3_hat']
    # if reference mu2 hat and reference mu3 hat is provided, use that instead.


    gwash_res = sims.do_GWASH_from_ldscores() # in practice, compute from ldscores only.
    if my_data_gen.ref_mu2_hat is not None:
        gwash_mu2_hat = my_data_gen.ref_mu2_hat

    if my_data_gen.ref_mu3_hat is not None:
        gwash_mu3_hat = my_data_gen.ref_mu3_hat

    GWASH_se_sample_theoretical = calc_gwash_se(my_data_gen.n,my_data_gen.m,gwash_mu2_hat,gwash_mu3_hat,gwash_res['gwash'])
    gwash_res_df = pd.DataFrame({"h2_gwash":gwash_res['gwash'],"h2_gwash_se": GWASH_se_sample_theoretical},index = [i])

    theory_res_df = pd.DataFrame({"h2_theory":h2_pop + (sigma_s**2 / Cf)},index = [i])

    res_dfs.append(pd.concat([ldsc_reg_res,ldsc_fixed_res,ldsc_fixed_nw_res,ldsc_nw_res,ldsc_nts_res,ldsc_nwts_res,gwash_res_df,theory_res_df],axis = 1))
res_dfs = pd.concat(res_dfs,axis = 0)
res_dfs.to_csv(DATA_DIR / "python_refpanel_indepf_n2000" / 'AR1_n{n}_m{m}_Fst{Fst_label}_sigmas03_rho03_h2pop03.csv'.format(n = n, m = m, Fst_label = Fst_label),index = None)
