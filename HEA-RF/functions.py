import numpy as np
import pandas as pd

def get_chemicals_list(filename):
    chemicals = pd.read_csv(filename)
    chemicals_list = list(chemicals[' symbol'].str.strip())
    return chemicals_list

def parse_chemical_name(name, chemicals_list):
    temp = ''
    stoi = []
    chem = []

    def cf(m):
        ls = []
        for sf in m:
            if sf.isupper():
                m = m.replace(sf," {}".format(sf))
        return m
    a = list(cf(name).split())

    for ele in a:
        if len(ele) == 1:
            stoi.append(1)
        elif ele[1].isdigit():
            stoi.append(float(ele[1:]))
        elif len(ele) > 2:
            stoi.append(float(ele[2:]))
        else:
            stoi.append(1)

        chem.append(''.join(i for i in ele if not (i.isdigit() or i ==".")))
    return chem, stoi

def check_artifacts_name(name, chem_list):
    if '(' and ')' in name:
        return filter_name(name, chem_list)
    else:
        return name

def filter_name(name, chem_list):
    output = ''
    temp = ''
    mult = 1
    i = 0
    j = 0
    count = 0
    k = 0
    for char in name:
        if char == '(':
            i = count
        if char == ')':
            j = count
            if j < len(name) - 1:
                if name[j+1].isdigit():
                    mult = name[j+1]
                    k = j + 1
        count += 1

    output = name[:i] + name[k+1:]

    temp = name[i+1:j]

    fix, stoi = parse_chemical_name(temp, chem_list)

    for el in fix:
        if int(mult) > 1:
            output += el + mult
        else:
            output += el

    return output


