import functions as f
import numpy as np
import math
import pandas as pd

class Features(object):

    def __init__(self, name, chem_list, data):
        self.name = name
        self.chem_list = chem_list
        self.data = data
        self.composition = np.zeros((len(self.chem_list), 1))
        self.weights = []
        self.ele = []
        self.stoi = []
        self.avg_en = 0
        self.avg_vec = 0
        self.avg_fbi = 0
        self.avg_xx = 0
        self.avg_elasticity = 0
        self.avg_space_group = 0
        self.avg_r = 0
        self.avg_r1 = 0
        self.avg_phi = 0
        self.avg_h = 0
        self.avg_meltingPoint = 0
        self.avg_ss = 0

    def elemental_composition(self):
        self.name = f.check_artifacts_name(self.name, self.chem_list)
        self.ele, self.stoi = f.parse_chemical_name(self.name, self.chem_list)

        for i in range(len(self.chem_list)):
            for j in range(len(self.ele)):
                if self.chem_list[i] == self.ele[j]:
                    self.composition[i, 0] = self.stoi[j]

        self.composition = self.composition/(np.sum(self.composition))
        self.weights = self.stoi/(np.sum(self.stoi))

        return self.composition, self.ele, self.stoi

    def get_rest_features(self):
        atomic_electroneg = self.data['electronegativity']

        atomic_fbi = self.data['fbi']
        atomic_space_group = self.data['space_group']
        atomic_r = self.data['radius']
        atomic_phi = self.data['Phi']
        atomic_elasticity = self.data['elasticity']
        atomic_meltingPoint = self.data['meltingPoint']

        s_electrons = self.data['s'].fillna(0.0)
        p_electrons = self.data['p'].fillna(0.0)
        d_electrons = self.data['d'].fillna(0.0)
        f_electrons = self.data['f'].fillna(0.0)

        for i in range(len(self.ele)):
            element = self.ele[i]
            self.avg_en += atomic_electroneg.loc[element]*self.weights[i]

            self.avg_vec += (s_electrons.loc[element]+p_electrons.loc[element]+
                             d_electrons.loc[element]+f_electrons.loc[element])*self.weights[i]
            self.avg_fbi += atomic_fbi.loc[element]*self.weights[i]
            self.avg_space_group += atomic_space_group.loc[element]*self.weights[i]
            self.avg_r += atomic_r.loc[element]*self.weights[i]
            self.avg_phi += atomic_phi.loc[element]*self.weights[i]
            self.avg_elasticity += atomic_elasticity.loc[element]*self.weights[i]
            self.avg_meltingPoint += atomic_meltingPoint.loc[element]*self.weights[i]
            self.avg_ss += math.log(self.weights[i])*self.weights[i]

        for i in range(len(self.ele)):
            element = self.ele[i]
            self.avg_xx += (atomic_electroneg.loc[element] - self.avg_en)**2*self.weights[i]
            self.avg_r1 += (1 - atomic_r.loc[element]/self.avg_r)**2*self.weights[i]

        hh  = pd.read_excel('H2.xlsx', index_col="element")
        for i in range(len(self.ele)-1):
            element_1 = self.ele[i]
            for j in range(i+1,len(self.ele)):
                element_2 = self.ele[j]
                self.avg_h += 4*hh.loc[element_1,element_2]*self.weights[i]*self.weights[j]

        self.avg_xx = self.avg_xx**(0.5)
        self.avg_r1 = self.avg_r1**(0.5)*100
        self.avg_ss = self.avg_ss*(-8.314)

        output = np.zeros(10)

        output[0] = self.avg_vec            ### VEC
        output[1] = self.avg_fbi            ### FBI
        output[2] = self.avg_xx             ### electronegativity
        output[3] = self.avg_elasticity     ### elasticity
        output[4] = self.avg_phi            ### work function
        output[5] = self.avg_r1             ### radius
        output[6] = self.avg_space_group    ### space group number
        output[7] = self.avg_h              ### enthalpy
        output[8] = self.avg_meltingPoint   ### meltingpoint
        output[9] = self.avg_ss             ### entropy

        return output