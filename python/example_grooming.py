import centauro 
import numpy as np
import math
import re
import os
import glob
import sys
import subprocess

# =============================================================================
#   Loading event information (particle momenta) ~ (p0,px,py,pz)
# =============================================================================
output_folder = "../output/example_grooming/"
mkdirCmd = ["mkdir",  output_folder]
subprocess.run(mkdirCmd)

z_cut   = 0.2  # grooming parameter

in_file = open("../input/EIC_eCM@63_test_12.txt" , "r")

event = 0
line = in_file.readline()
clean_line = np.array(re.sub(' +', ' ', line.strip()).split(" "))
while (clean_line[0] == "<event_init>"): 
    Q = clean_line[1].astype(np.float)
    x = clean_line[2].astype(np.float)
    event += 1
    particles = []
    line = in_file.readline()
    clean_line = np.array(re.sub(' +', ' ', line.strip()).split(" "))
    while (clean_line[0] != "<event_end>"):
        aline = clean_line.astype(np.float)
        particles.append(aline.tolist())
        line = in_file.readline()
        clean_line = np.array(re.sub(' +', ' ', line.strip()).split(" "))

    line = in_file.readline()
    clean_line = np.array(re.sub(' +', ' ', line.strip()).split(" "))    
    n_particles = len(particles)    
    
    if (n_particles <= 0): 
        print(n_particles)
        continue


    jets_init = []
    inclusive_jets = []
    

    rec = 1
    for particle in particles:
        jets_init.append(centauro.jet([particle], rec, rec ))
        rec += 1

    out_file = open(output_folder + "ungroomed_eCM@63_{0}_HM.txt".format(event), "w")
    np.savetxt(out_file, np.array(particles) , delimiter = "\t")
    out_file.close()

    groomed = centauro.centauro_grooming(jets_init, z_cut)
    
    out_file = open(output_folder+"groomedCentauro_eCM@63_{0}_HM.txt".format(event), "w")
    np.savetxt(out_file, np.array(groomed.constituents) , delimiter = "\t")
    out_file.close()

    if (event % 1000 == 0): print("{0} events loaded and analyzed.".format(event))



