import os
import multiprocessing

import h5py
import numpy as np
import pandas as pd

def distancefromcentre(cx, cy, cz, x, y, z, ):
    """
    

    Parameters
    ----------
    cx : Halo centre x-coord float
    cy : Halo centre y-coord float
    cz : Halo centre z-coord float
    x : Particle Position x-coord
    y : Particle Position y-coord
    z : Particle Position z-coord

    Returns
    -------
    Distance between particle and halo centre

    """
    return (np.sqrt((np.power(np.subtract(x,cx), 2)+ np.power(np.subtract(y,cy), 2) + np.power(np.subtract(z,cz), 2)))) # distance between the centre and given point





n_process = 50
box_size = 50
run = 1

lhalotree_dir = f'/disk01/rmcg/downloaded/tng/tng{box_size}-{run}/merger_tree/lhalotree/'

def extract_trees(filepath):
    print(f'Starting processing file: {filepath}')
    with h5py.File(filepath, 'r') as file:
        n_halos_in_tree = np.array(file['/Header/TreeNHalos'])

        trees = []
        subhalo_ids = []
        for i_tree, n_halo in enumerate(n_halos_in_tree):
            arr = {}
            tree = file[f'Tree{i_tree}']

            # TODO: Convert mass units
            arr_mass_type = np.array(tree['SubhaloMassType'])
            arr_position = np.array(tree['SubhaloPos'])
            arr['positionX'] = arr_position[:, 0]
            arr['positionY'] = arr_position[:, 1]
            arr['positionZ'] = arr_position[:, 2]
            arr['halfmass_rad'] = np.array(tree['SubhaloHalfmassRad'])
            arr['mass'] = np.array(tree['SubhaloMass'])
            arr_spin = np.array(tree['SubhaloSpin'])
            arr['spinX'] = arr_spin[:, 0]
            arr['spinY'] = arr_spin[:, 1]
            arr['spinZ'] = arr_spin[:, 2]
            arr['bh_mass'] = np.array(tree['SubhaloBHMass'])
            arr['bh_dot'] = np.array(tree['SubhaloBHMdot'])
            arr['dm_mass'] = arr_mass_type[:, 1]
            arr['gas_mass'] = arr_mass_type[:, 0]
            #arr['gas_metallicity'] = np.array(tree['SubhaloGasMetallicity'])
            arr['sfr'] = np.array(tree['SubhaloSFR'])
            arr['stellar_mass'] = arr_mass_type[:, 4]
            arr['vel_dispersion'] = np.array(tree['SubhaloVelDisp'])
            arr['v_max'] = np.array(tree['SubhaloVmax'])
            #arr['stellar_metallicity'] = np.array(tree['SubhaloStarMetallicity'])
            arr['particle_number'] = np.array(tree['SubhaloLen'])
            
            arr['main_prog_index'] = np.array(tree['FirstProgenitor'])
            arr['snap_num'] = np.array(tree['SnapNum'])
            arr['subhalo_id'] = np.array(tree['SubhaloNumber'])
            
            fof_number = np.array(tree['SubhaloGrNr'])
            radial_distance = distancefromcentre(tree['GroupPos'][fof_number][:, 0], 
                                                 tree['GroupPos'][fof_number][:, 1], 
                                                 tree['GroupPos'][fof_number][:, 2], 
                                                 arr_position[:, 0], arr_position[:, 1], 
                                                 arr_position[:, 2])
            arr['fof_mass'] = np.array(tree['GroupMass'][fof_number])
            arr['fof_distance'] = radial_distance
            
            arr_central_index = np.array(tree['FirstHaloInFOFGroup'])
            arr['is_central'] = np.zeros(n_halo, dtype=bool)
            for i_halo, i_central in enumerate(arr_central_index):
                arr['is_central'][i_halo] = (i_halo == i_central)
            
            min_snap = 2
            max_snap = 99
            input_properties = ['positionX','positionY','positionZ','halfmass_rad','mass',
                                'particle_number','gas_mass','dm_mass','stellar_mass', 'bh_mass', 
                                'spinX','spinY','spinZ', 'vel_dispersion','v_max', 'bh_dot',
                                'sfr', 'fof_mass','fof_distance']

            snapshots = list(range(max_snap, min_snap-1, -1))
            n_input, n_snap = len(input_properties), len(snapshots)
            input_features = [str(snap)+prop for snap in snapshots for prop in input_properties]

            valid = (arr['snap_num'] == max_snap)
            # TODO: Set mass cut
            valid &= (arr['particle_number'] > 500)
            n_valid_sub_this_file = np.sum(valid)
            if n_valid_sub_this_file == 0:
                continue

            i_sub = 0
            histories = np.zeros((n_valid_sub_this_file, n_input*n_snap), dtype='float64')
            for i_halo in np.arange(n_halo)[valid]:
                i_prog = i_halo
                while i_prog != -1:
                    snap_num = arr['snap_num'][i_prog]
                    if snap_num < min_snap:
                        break
                    
                    positionX = arr['positionX'][i_prog]
                    positionY = arr['positionY'][i_prog]
                    positionZ = arr['positionZ'][i_prog]
                    halfmass_rad = arr['halfmass_rad'][i_prog]
                    mass = arr['mass'][i_prog]
                    particle_number = arr['particle_number'][i_prog]
                    bh_mass = arr['bh_mass'][i_prog]
                    bh_dot = arr['bh_dot'][i_prog]
                    dm_mass = arr['dm_mass'][i_prog]
                    gas_mass = arr['gas_mass'][i_prog]
                    spinX = arr['spinX'][i_prog]
                    spinY = arr['spinY'][i_prog]
                    spinZ = arr['spinZ'][i_prog]
                    vel_dispersion = arr['vel_dispersion'][i_prog]
                    v_max = arr['v_max'][i_prog]
                    fof_mass = arr['fof_mass'][i_prog]
                    fof_distance = arr['fof_distance'][i_prog]
                    #gas_metallicity = arr['gas_metallicity'][i_prog]
                    sfr = arr['sfr'][i_prog]
                    stellar_mass = arr['stellar_mass'][i_prog]
                    #stellar_metallicity = arr['stellar_metallicity'][i_prog]

                    i_start = (max_snap - snap_num) * n_input
                    # This has to line up with where input columns are defined
                    #data = [bh_mass, bh_dot, dm_mass, gas_mass, gas_metallicity,
                    #        sfr, stellar_mass, stellar_metallicity]
                    data = [positionX,positionY,positionZ,halfmass_rad,mass,
                                        particle_number,gas_mass,dm_mass,stellar_mass, bh_mass, 
                                        spinX,spinY,spinZ, vel_dispersion,v_max, bh_dot,
                                        sfr, fof_mass,fof_distance]
                    histories[i_sub, i_start:i_start+n_input] = data

                    i_prog = arr['main_prog_index'][i_prog]

                i_sub += 1

            trees.append(pd.DataFrame(histories, columns=input_features))
            subhalo_ids.append(arr['subhalo_id'][i_halo])
    return pd.concat(trees, ignore_index=True), subhalo_ids


pool_result = []
pool = multiprocessing.Pool(n_process)
filenames = [lhalotree_dir+name for name in os.listdir(lhalotree_dir)]
while filenames:
    files_to_process, filenames = filenames[:n_process], filenames[n_process:]
    pool_result += pool.map(extract_trees, files_to_process)

all_trees, all_subhalo_ids = pool_result.pop(0)
while pool_result:
    trees, subhalo_ids = pool_result.pop(0)
    all_trees = pd.concat([all_trees, trees], ignore_index=True)
    all_subhalo_ids += subhalo_ids

# TODO: Save data
print(all_trees)
print(all_subhalo_ids)
all_trees.to_csv('50-1-subhalo-history.csv',index_label=all_subhalo_ids)
