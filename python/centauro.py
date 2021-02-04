import numpy as np
import math
import re
import os
import glob
import sys

# Four-vector functons 
# --------------------

def length(vec4): return sum([x**2 for x in vec4[1:4]])**0.5

def pT(vec4):     return sum([x**2 for x in vec4[1:3]])**0.5

def phi(vec4):
    e = 0.000000000001
    return math.copysign(math.acos(vec4[1]/((vec4[1]**2 + vec4[2]**2)**0.5 + e)), vec4[2])

def eta(vec4): 
    e = 0.000000000001
    return math.atanh( vec4[3]/(length(vec4)+e) )

def theta(vec4_1, vec4_2 = np.array([1,0,0,1])): 
    e = 0.000000000001
    vec3_1 = vec4_1[1:4]
    vec3_2 = vec4_2[1:4]
    return  math.acos( np.dot(vec3_1, vec3_2) / ( (length(vec4_1)+e) * (length(vec4_2)+e) ) )

def dot(vec4_1, vec4_2):
    return vec4_1[0] * vec4_2[0] - (vec4_1[1] * vec4_2[1]+ vec4_1[2] * vec4_2[2]+vec4_1[3] * vec4_2[3])

def mass(vec4): 
    m2 = vec4[0]**2 - (vec4[1]**2 + vec4[2]**2 +vec4[3]**2 ) 
    if (m2 >= 0): return  math.sqrt(m2)
    else: return 0

# Class for pseudo-jets 
# ---------------------
    
class jet:
    
    def __init__(self,  constituents, record, daughter_1, daughter_2 = 0): 
        self.constituents = constituents
        self.record = record 
        self.daughter_1 = daughter_1 
        self.daughter_2 = daughter_2  
        
    #-----------------------------------  
    # n    = [1, 0, 0, +1] : pP = p.n
    # nbar = [1, 0, 0, -1] : pM = p.nbar
    #-----------------------------------
    # static methods
    @staticmethod
    def total_p(partices):return [sum(x) for x in zip(*partices)]
    @staticmethod
    def length(four_vector): return sum([x**2 for x in four_vector[1:4]])**0.5
    @staticmethod
    def pL(four_vector): return four_vector[3]  
    @staticmethod
    def pM_SM(four_vector): return four_vector[0] + four_vector[3]
    @staticmethod
    def pP_SM(four_vector): return four_vector[0] - four_vector[3]
    #-----------------------------------
    def p(self):
        return self.total_p(self.constituents)

    def p3(self):
        vec4 = self.total_p(self.constituents)
        return vec4[1:4]
    
    def e(self): 
        four_vector = self.total_p(self.constituents)
        return four_vector[0]
    
    def pT(self): 
        four_vector = self.total_p(self.constituents)
        return sum([x**2 for x in four_vector[1:3]])**0.5
    
    def pM(self): 
        four_vector = self.total_p(self.constituents)
        return self.pM_SM(four_vector)
    
    def pP(self):
        four_vector = self.total_p(self.constituents)
        return self.pP_SM(four_vector)
    
    def eta(self): 
        four_vector = self.total_p(self.constituents)
        e = 0.000000000001
        return math.atanh(self.pL(four_vector)/(self.length(four_vector)+e))
    
    def etaB(self): 
        four_vector = self.total_p(self.constituents)
        e = 0.000000000001
        pT = sum([x**2 for x in four_vector[1:3]])**0.5
        return 2 * pT / (self.pP_SM(four_vector) + e )

    def etaBP(self): 
        four_vector = self.total_p(self.constituents)
        e = 0.000000000001
        pT = sum([x**2 for x in four_vector[1:3]])**0.5
        return math.asinh( 2 * pT / (self.pP_SM(four_vector) + e ) )
    
    def phi(self):
        four_vector = self.total_p(self.constituents)
        e = 0.000000000001
        return math.atan(four_vector[2]/(four_vector[1] + e))
        
    def y(self): 
        four_vector = self.total_p(self.constituents)
        e = 0.000000000001
        return math.log(self.pM_SM(four_vector)/(self.pP_SM(four_vector)+e))/2

    def mass2(self):
        vec4 = self.p()
        return dot( vec4, vec4 )
    #-----------------------------------
    
    def join(self,other,rec): 
        return jet(self.constituents + other.constituents, rec, self.record, other.record )
    #-----------------------------------
    
    def __add__(self, other):
        return jet(self.constituents + other.constituents, 0, 0)
    #-----------------------------------
        
    def __str__(self):
        return "pseudo-jet record number:{0}  daughters:({1},{2})".format(self.record, self.daughter_1, self.daughter_2)

    def thrust_A(self, Q):
        thrustTot = 0
        particles = self.constituents
        vec4 = self.p()
        norm = length(vec4)
        axis = [1, vec4[1]/norm, vec4[2]/norm, vec4[3]/norm ]
        for particle in particles:
            thrustTot += min(dot(particle, axis), self.pP_SM(particle))
        return thrustTot / Q 

    def thrust_B(self, Q):
        thrustTot = 0
        particles = self.constituents
        for particle in particles:
            thrustTot += min(self.pM_SM(particle), self.pP_SM(particle))
        return thrustTot / Q 
        
    #-----------------------------------
    #           END OF CLASS
    #-----------------------------------
    

pi = math.pi 
e  = math.e


#------------------------------------------------------------------------------   
def DR2(jet1, jet2):
    eta1 = jet1.etaBP()
    eta2 = jet2.etaBP()
    phi1 = jet1.phi()
    phi2 = jet2.phi()
    return  ( (eta1 - eta2)**2 + 2 * eta1 * eta2 *(1- math.cos(phi1-phi2) ) ) 
#------------------------------------------------------------------------------
def d_init(jets ):
    d = []
    jet_N = len(jets)
    for i in range(0,jet_N):
        row = []
        for j in range(0, jet_N):
            dr2 = DR2(  jets[i] , jets[j] )
            row.append(dr2 )
        d.append(row)
    return d
#------------------------------------------------------------------------------
def d_next(jets, loc ):
    jet_N = len(jets)
    row = []
    for j in range(0,jet_N):
        dr2 = DR2(  jets[loc] , jets[j] )
        row.append( dr2 )
    return row
#------------------------------------------------------------------------------
def dF_init(jets ):
    d_flatten = []
    jet_N = len(jets)
    for i in range(0,jet_N-1):
        for j in range(i+1,jet_N):
            dr2 = DR2(  jets[i] , jets[j] )
            d_flatten.append(dr2)
    return d_flatten
#------------------------------------------------------------------------------
def find_index(N_jets, Y):
    M = 0
    N = N_jets - 1
    while ( (Y -N + M) >= 0 ): 
        Y = Y - N + M
        M += 1
    return [M, M + 1 + Y]
#------------------------------------------------------------------------------
def flatten(list):
    a = []
    jet_N = len(list)
    for i in range(0,jet_N-1):
        for j in range(i+1,jet_N):
            a.append(list[i][j])
    return a
#------------------------------------------------------------------------------
    
def centauro_clustering(jets_in, R = 1):
    
    jets = jets_in[:] 
    inclusive_jets = []
    
    d_list  = d_init(jets)
    dF_list = dF_init(jets)
    
    while (len(jets) > 1 ):
        dF_min = min( dF_list )
        if ( dF_min > R**2 ):
            n_jets = len(jets)
            for loc in range(0, n_jets):
                inclusive_jets.append( jets[loc] )
            return (inclusive_jets)

        else:
            n_jets = len(d_list)
            Y = dF_list.index(dF_min)
            loc = find_index(n_jets, Y)
            jets[loc[0]] =  jets[loc[0]] + jets[loc[1]]
            del jets[loc[1]]
            del d_list[loc[1]]
            for i in range(0,len(d_list)):
                del d_list[i][loc[1]]
            new_row = d_next(jets,loc[0])
            d_list[loc[0]] = new_row
            for i in range(0,len(d_list)):
                d_list[i][loc[0]] = new_row[i]
        dF_list = flatten(d_list)
        
    inclusive_jets.append( jets[0] )
    
    return inclusive_jets


#------------------------------------------------------------------------------
    
def clustering_No_beam(jets_in):
    

    jets = jets_in[:] 
    inclusive_jets = []

    rec = len(jets_in) + 1

    first_beam_jet = True
    
    d_list  = d_init(jets)
    dF_list = dF_init(jets)

    if (len(jets) > 1 ):
    
        while (len(jets) > 1 ):
            dF_min = min( dF_list )

            n_jets = len(d_list)
            Y = dF_list.index(dF_min)
            loc = find_index(n_jets, Y)
            jets[loc[0]] =  (jets[loc[0]]).join( jets[loc[1]], rec  )
            rec += 1
            inclusive_jets.append( jets[loc[0]] )
            del jets[loc[1]]
            del d_list[loc[1]]
            for i in range(0,len(d_list)):
                del d_list[i][loc[1]]
            new_row = d_next(jets,loc[0])
            d_list[loc[0]] = new_row
            for i in range(0,len(d_list)):
                d_list[i][loc[0]] = new_row[i]
            dF_list = flatten(d_list)
    else:
        inclusive_jets.append( jets[0] )
    
    return inclusive_jets


def De_clustering(tree,  z_cut):

    fail  = True 
    tree_index = -1

    while(fail == True):
        d1 = tree[tree_index].daughter_1 
        d2 = tree[tree_index].daughter_2 
        if (d1 == 0 or d2 == 0 ): break
        z1 =  tree[d1-1].pP() 
        z2 =  tree[d2-1].pP() 
        zT = z1 + z2
        r = min(z1, z2)/zT
        if (r < z_cut) : 
            if (z1 < z2 ): tree_index = d2 -1
            else: tree_index = d1 -1
        else: break
    return tree[tree_index]

def centauro_grooming(event, z_cut):

    jets_init = event

    jets_temp = event

    clust =  jets_init + clustering_No_beam(jets_temp) 

    return De_clustering(clust,  z_cut)

#------------------------------------------------------------------------------