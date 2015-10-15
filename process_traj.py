#!/bin/env python
import numpy as np
import mdtraj as mdt
from mdtraj import element
from joblib import Parallel, delayed  
import multiprocessing
import os

import glob
import mdtraj as md
from msmbuilder.utils import verbosedump,verboseload
from msmbuilder.featurizer import DihedralFeaturizer
import wetmsm as wmsm

#load the trajset

def calc_obs(traj):
    arg_cz_id = 2442
    glu_cd_id = 862
    lys_nz_id = 634
    tyr_oh_id = 2019
    inactive = mdt.load("./topologies/inactive.pdb")
    active = mdt.load("./topologies/active.pdb")    
    
    aloop_atoms_list=[i.index for residue in np.arange(147,168) for i in inactive.topology.residue(residue).atoms]
    all_heavy = [i.index for i in inactive.topology.atoms if i.residue.is_protein and
 i.element.name != "hydrogen"]
    print("Processing %s"%traj)
    #load the trajectory
    trj = mdt.load(traj,atom_indices=np.arange(inactive.n_atoms))

    inactive_rms = mdt.rmsd(trj,inactive,atom_indices=all_heavy)
    active_rms = mdt.rmsd(trj,active,atom_indices = all_heavy)
    aloop_rms = mdt.rmsd(trj,inactive,frame=0,atom_indices=aloop_atoms_list)
    distances= mdt.compute_distances(trj,np.vstack(([arg_cz_id,glu_cd_id],       [lys_nz_id,glu_cd_id])))
    return dict(fname=os.path.basename(traj),inactive_rmsd=inactive_rms,active_rmsd=active_rms,aloop_inactive_rmsd=aloop_rms,glu_arg = distances[:,0],gly_lys=distances[:,1])



def feat_traj(traj):
    #load again to get the waters
    trj = mdt.load(traj)
    atp_solute=[i.index for i in trj.topology.atoms if (i.residue.name == "atp" and (i.element.name == "oxygen" or i.element.name == "nitrogen")) or (i.residue.name == "MG" )]
     #get the oxygen and nitrogen indices
    solute_indices = [i.index for i in trj.topology.atoms if i.residue.is_protein and (i.element.name=="oxygen" or i.element.name=="nitrogen")]
     #get the oxygen solvent indices 
    solvent_indices = [i.index for i in trj.topology.atoms if (i.residue.is_water and i.element.name != "hydrogen")]
     
    #set up featurizers 
    atp_feat = wmsm.SolventShellsFeaturizer(atp_solute,solvent_indices,2,0.3)
    water_feat = wmsm.SolventShellsFeaturizer(solute_indices,solvent_indices,2,0.3)
    dihedral_feat = DihedralFeaturizer(["phi","psi","chi1"])
    
    #calculate features
    water_features = water_feat.partial_transform(trj)
    dihedral_features = dihedral_feat.partial_transform(trj)
    atp_features = atp_feat.partial_transform(trj)
      
    combined_features = np.hstack((dihedral_features,water_features,atp_features))
     
    return combined_features
     #dump
    
    fname = os.path.basename(traj)
    save_path = os.path.join("/nobackup/msultan/research/kinase/fyn_kinase/fah_data/features/")
    verbosedump(dihedral_features,os.path.join((save_path,"dihedral/%s"%fname)))
    verbosedump(water_features,os.path.join((save_path,"water/%s"%fname)))
    verbosedump(atp_features,os.path.join((save_path,"atp/%s"%fname)))
    verbosedump(combined_features,os.path.join((save_path,"combined/%s"%fname)))
    return 


#verbosedump(result_dict,'scatter.h5')
