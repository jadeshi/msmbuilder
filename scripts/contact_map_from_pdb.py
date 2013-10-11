#!/usr/bin/env python
from msmbuilder import Trajectory
import numpy as np
from msmbuilder.geometry import contact
from msmbuilder import arglib


def get_residue_contacts(traj):
    '''Returns a 2D numpy array enumerating all unique pairs of residues given a PDB loaded in trajectory format.'''
    residues = np.unique(traj['ResidueID'])
    residue_contacts = np.array([])
    for i in range(len(residues)):
        for j in range(len(residues)):
            if j >= i and j != i:
                pair = np.array([i, j])
                if len(residue_contacts) == 0:
                    residue_contacts = pair
                else:
                    residue_contacts = np.vstack((residue_contacts, pair))
    return residue_contacts


def get_residue_membership(traj):
    '''Get the "residue membership" list of lists that is required for running contact.py.'''
    resids = traj['ResidueID']
    residues = np.unique(resids)
    residue_membership = []
    for i in residues:
        residue_membership.append(np.where(resids == i)[0])
    return residue_membership


def get_contact_distances(traj, residue_membership, residue_contacts):
    '''Uses contact.py to get the contact distance between all pairs of residues specified by "residue_contacts", and returns a n_contacts x 3 numpy array with the first 2 columns specifying the residue pair and the 3rd the contact distan
ce, in nanometers.'''
    xyzlist = traj['XYZList']
    contact_distances = np.array([contact.residue_distances(xyzlist, residue_membership, residue_contacts)[0]])
    return np.concatenate((residue_contacts.transpose(), contact_distances)).transpose()


def get_contact_map(contact_distances, cutoff):
    '''Converts contact_distances (see "get_contact_distances") into a contact map (n_contacts x 3 numpy array).'''
    for i in range(len(contact_distances)):
        if contact_distances[i][2] < cutoff:
            contact_distances[i][2] = 1
        else:
            contact_distances[i][2] = 0
    return contact_distances


def contact_map_to_matrix(contact_map, read=True):
    '''Converts a contact map into a square matrix. Can either be a 3-column text file (default) or a n_contacts x 3 numpy array. If the "read" argument is set to not "True", then contact map is assumed to be a numpy array.'''
    if read == True:
        contact_map = np.loadtxt(contact_map)
    n = (1 + (1 + 4 * 2 * len(contact_map)) ** 0.5) / 2
    matrix = np.zeros((n, n))    
    for contact in contact_map:
        matrix[contact[0]][contact[1]] = contact[2]
    return matrix
        

if __name__ == '__main__':
    parser = arglib.ArgumentParser(description='''Creates a contact map given a PDB and a contact cutoff. The contact map will be a 3-column text file, with the first two columns listing the residue pair, and the third column being either
 0 if the two residues are not in contact or 1 if they are. The map will enumerate all unique residue pairs once, and excludes self-contacts. The text map can be read and converted into a numpy square matrix format, with the [i,j]th entry
 corresponding to contact between residues i and j, using the utility function "contact_map_to_matrix". If desired, the contact distances between residue pairs can be explicitly written to a separate 3-column text file, with the first two
 columns listing the residue pair and the 3rd giving the contact distance in nanometers.''')
    parser.add_argument('input_PDB', type=str, help='''The PDB from which you want to extract a contact map.''')
    parser.add_argument('contact_cutoff', type=str, help='''Cutoff you want to use to define a contact, in nanometers.''', default='0.3')
    parser.add_argument('output_map', type=str, help='''The name you want to give the contact map text file.''')
    parser.add_argument('output_distances', type=str, help='''The name you want to give the contact distances text file. This file will not be written unless an argument is supplied.''', default='-1')
    args = parser.parse_args()
    
    PDB = Trajectory.load_trajectory_file(args.input_PDB)
    residue_membership = get_residue_membership(PDB)
    residue_contacts = get_residue_contacts(PDB)
    contact_distances = get_contact_distances(PDB, residue_membership, residue_contacts)
    if args.output_distances != '-1':
        np.savetxt(args.output_distances, contact_distances, fmt='%i %i %f')
    contact_map = get_contact_map(contact_distances, float(args.contact_cutoff))
    np.savetxt(args.output_map, contact_map, fmt='%i %i %i')
