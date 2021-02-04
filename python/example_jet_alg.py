import numpy as np
import math
import re
import os
import glob
import sys
import centauro 
import subprocess

# =============================================================================
#   Loading event information (particle momenta) ~ (p0,px,py,pz)
# =============================================================================

output_folder = "../output/example_jet_alg/"
mkdirCmd = ["mkdir",  output_folder]
subprocess.run(mkdirCmd)


R = 0.8

in_file  = open("../input/single_event_example.txt","r")

event = 0
line = in_file.readline()
clean_line = np.array(re.sub(' +', ' ', line.strip()).split(" "))
while (clean_line[0] == "<event_init>"): 
    event += 1
    particles = []
    line = in_file.readline()
    clean_line = np.array(re.sub(' +', ' ', line.strip()).split(" "))
    while (clean_line[0] != "<event_end>"):
        aline = clean_line.astype(np.float)
        particles.append(aline.tolist())
        line = in_file.readline()
        clean_line = np.array(re.sub(' +', ' ', line.strip()).split(" "))
        
    n_particles = len(particles)    
    

    jets = []
    inclusive_jets = []
    
    rec = 1
    for particle in particles:
        jets.append(centauro.jet([particle], rec, rec))
            
    inclusive_jets = centauro.centauro_clustering(jets, R)

    i = 0
    for jet_i in inclusive_jets:
        out_file = open(output_folder + "jet_{0}.txt".format(  i ) ,"w")
        for cont_i in  jet_i.constituents:
            vector = [[ centauro.eta(cont_i), centauro.phi(cont_i), centauro.theta(cont_i), cont_i[0] ]]
            np.savetxt(out_file, vector , delimiter = "\t")
        i += 1
        out_file.close()
    

    line = in_file.readline()
    clean_line = np.array(re.sub(' +', ' ', line.strip()).split(" "))

print("{0} events loaded and analyzed.".format(event))
    
out_file.close()

