#!/bin/env/python

import mdtraj as mdt
import os
from IPython import parallel
import numpy as np
from fah_reseeder.convert_project import *
import wetmsm as wmsm
from msmbuilder.utils import verboseload,verbosedump
from msmbuilder.decomposition import tICA
import glob
import pandas as pd

def keynat(string):
    '''A natural sort helper function for sort() and sorted()
    without using regular expression.
    >>> items = ('Z', 'a', '10', '1', '9')
    >>> sorted(items)
    ['1', '10', '9', 'Z', 'a']
    >>> sorted(items, key=keynat)
    ['1', '9', '10', 'Z', 'a']
    '''
    r = []
    for c in string:
        try:
            c = int(c)
            try:
                r[-1] = r[-1] * 10 + c
            except:
                r.append(c)
        except:
            r.append(c)
    return r

def three_to_one(amino):
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    return d[amino]

def featurize_traj(job_tuple):
    #separate out the job tuple into required things
    mutant,mutant_dir,project,proj_folder,proj_top_folder,traj_file,stride,save_common,allowed_residue_ind \
    = job_tuple
    #load top file to setup solute/solvent indices
    top_path = os.path.join(proj_top_folder, "%s.pdb"%os.path.basename(traj_file).split("_")[0])
    top_trj = mdtraj.load(top_path)

    #set up featurizer objects
    dihedral_feat = DihedralFeaturizer(types=['phi', 'psi','chi1'])

    #load the trajectory
    try:
        trj = mdtraj.load(traj_file,stride=stride)
    except:
        print "Cant featurize %s"%traj_file
        return 
    #setup file name
    traj_name = os.path.splitext(os.path.basename(traj_file))[0]
    print traj_name
    dihedral_output_file = os.path.join(mutant_dir,"features/dihedral_features/")+str(project)+\
    "_"+traj_name+".h5"
    water_output_file = os.path.join(mutant_dir,"features/water_features/")+str(project)+\
    "_"+traj_name+".h5"
    combined_output_file = os.path.join(mutant_dir,"features/combined_features/")+str(project)+\
    "_"+traj_name+".h5"
    do_again=True
    already_done=False
    if os.path.isfile(combined_output_file):
    	f = verboseload(combined_output_file)
	if f.shape[0]!=trj.n_frames:
		already_done=True

    if not already_done or do_again:
        dihedral_features = dihedral_feat.partial_transform(trj)

    	traj_name = os.path.splitext(os.path.basename(traj_file))[0]

        dihedral_output_file = os.path.join(mutant_dir,"features/dihedral_features/")+str(project)+\
        "_"+traj_name+".h5"

        #now we can dump
    	verbosedump(dihedral_features,dihedral_output_file)

        if save_common:
            dih_df = pandas.DataFrame(dihedral_feat.describe_features(top_trj))

            dih_f_ind = numpy.array([set(i).issubset(allowed_residue_ind) for i in dih_df["resid"]])

            subset_dihedral_features = dihedral_features[:,dih_f_ind]

            dihedral_output_file = os.path.join(mutant_dir,"features/common_basis/dihedral_features/")+\
            str(project)+"_"+traj_name+".h5"


            #now we can dump
            verbosedump(subset_dihedral_features,dihedral_output_file)
            #save the featurizer information.
            verbosedump([dih_df,allowed_residue_ind,dih_f_ind,],\
os.path.join(mutant_dir,"features/common_basis/dihedral_features/")+"saved_dihed_feat.h5")

            return

    else:
	   print "skipping featurization for %s since its already done"%traj_name
    return


def featurize_project(mutant,project,proj_folder,proj_top_folder,mutant_dir,stride,save_common,allowed_residue_ind,view):
    feature_dict={}

    traj_list =  glob.glob(proj_folder+"/trajectories/*.hdf5")
    traj_list = sorted(traj_list, key=keynat)
    if save_common:
        #get common sequences
        #make jobs
        jobs = [(mutant,mutant_dir,project,proj_folder,proj_top_folder,traj_file,stride,\
            save_common,allowed_residue_ind ) for traj_file in traj_list]
    else:
        jobs = [(mutant,mutant_dir,project,proj_folder,proj_top_folder,traj_file,stride,None,None ) for \
    traj_file in traj_list]

    print "Featurizing %d trajectories for project %d"%(len(jobs),project)

    results = view.map_sync(featurize_traj,jobs)
        #print np.shape(result[0]),np.shape(result[1]),np.shape(np.hstack(result))
    return



def main():

    client_list = parallel.Client()
    print("Running on:",len(client_list.ids))
    view = client_list.direct_view()
    view.block = True
    with view.sync_imports():
        import numpy
        import mdtraj
        import wetmsm
        import os
        import pandas
        from msmbuilder.featurizer import DihedralFeaturizer
        from msmbuilder.utils import verbosedump,verboseload

    if os.getcwd()!="/nobackup/msultan/research/kinase/fyn_kinase/fah_data/PROJ9116":
        os.chdir("/nobackup/msultan/research/kinase/fyn_kinase/fah_data/PROJ9116")

    base_dir = os.getcwd()
    
    proj_folder = base_dir
    proj_top_folder = os.path.join((base_dir+"/topologies/"))
    print proj_folder
    print proj_top_folder
    extract_project_wrapper(proj_folder,proj_top_folder,view)
    
    #featurize_project(mutant,project,proj_folder,proj_top_folder,\
     #           mutant_dir,5,save_common,allowed_residue_ind,view)


if __name__ == '__main__':
    main()
