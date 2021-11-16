""""
This code uses the Monte Carlo Direct Method to simulate step polymerization
"""

# Imports --------------------------------------------------------------------------------------------------------------
import Paths
import random
import re
import tkinter as tk
import math
import itertools
import matplotlib.pyplot as plt
from molmass import Formula
from collections import Counter
import json
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import datetime
from typing import Union
import operator
import copy
import pickle
import os


# Global Variables -----------------------------------------------------------------------------------------------------

# Simulation
t_max = 10**30
t = 0

# Avogadro
Na = 6.03*10**23
R = 8.314    # J/mol.K

# Reading the parameters file
with open("Parameters_FromMonomer_Step.json", "r") as f:
    input_parameters = json.load(f)
    number_of_molecules = input_parameters["number_of_molecules"]
    initial_concentration = input_parameters["initial_concentration"]  # mol/m3
    molecules = input_parameters["initial_molecules"]
    T = input_parameters["temperature"]  # K
    kinetics = input_parameters["kinetics"]
    byproduct = input_parameters["byproduct"]

    # Reading types of endgroups
    smiles_a = molecules[0]["endgroup"]
    order_a = molecules[0]["order"]
    smiles_b = molecules[1]["endgroup"]
    order_b = molecules[1]["order"]

# Initial monomers
initial_monomers = [x["molecule"] for x in molecules]

# Output folder
output_folder = Paths.Directories("Step").output_folder

# Functions ------------------------------------------------------------------------------------------------------------
class EndSimulationException(Exception):
    # Ends the simulation
    pass

class NextIterationException(Exception):
    # Skips to next iteration
    pass

class TypeA:
    # Classifies start/end of molecule as A
    smiles = smiles_a
    order = order_a

class TypeB:
    # Classifies start/end of molecule as B
    smiles = smiles_b
    order = order_b

class ByProduct:
    # Classifies substance as byproduct
    smiles = byproduct
    mass = Formula(byproduct).mass

class Polymer:
    # Classifies substance as polymer
    pass

class Molecule:
    def __init__(self, smiles, start_type, end_type, mol_type=Polymer()):
        self.smiles = smiles
        self.start_type = start_type
        self.end_type = end_type
        self.type = mol_type

class System_FromMonomer:

    def __init__(self, entry_molecules_smiles: dict, molecule_type: dict, t: float, t_max: float, kinetics: dict, V=0.00001):
        self.t = t
        self.t_max = t_max
        self.initial_state = entry_molecules_smiles
        self.k = kinetics
        self.V = V
        self.state = self.initial_state
        self.rate_of_reaction_dict = None
        self.reaction_probability_dict = None
        self.pair_dict = None
        self.chosen_reaction = None
        self.tau = None
        self.reaction_type = None
        self.molecule_type = molecule_type
        # Time evolution
        self.time_evolution = {self.t: copy.deepcopy(self.state)}

    def rate_of_reaction(self, mol1: str, mol2: str) -> float:
        """
        Calculates the rate of reaction
        :param mol1: molecule 1
        :param mol2: molecule 2
        :return: rate of reaction
        """

        A = self.k["a"]    # L/mol.s
        Ea = self.k["Ea"]    # J/mol
        k = A * math.exp(-Ea / (R * T))
        if mol1 == mol2 or mol1 == mol2[::-1]:
            rate = k * 2 / (self.V * Na) * self.state[mol1] * (self.state[mol1] - 1) / 2
        else:
            rate = k / (self.V * Na) * self.state[mol1] * self.state[mol2]
        return rate

    def reverse_reaction(self, smiles1: str, smiles2: str) -> dict:
        """
        Generates the product of the reverse reaction
        :param smiles1: molecule 1
        :param smiles2: molecule 2
        :return: product
        """

        # Get the molecules
        mol1 = self.molecule_type[smiles1]
        mol2 = self.molecule_type[smiles2]

        if isinstance(mol1.type, ByProduct):
            _byprod = smiles1
            _mol = smiles2
        else:
            _byprod = smiles2
            _mol = smiles1


        # Choosing a product candidate
        _mols = [x for x in self.state.keys() if (len(x) < len(_mol))
                 and (not isinstance(self.molecule_type[x].type, ByProduct))
                 and (type(self.molecule_type[x].end_type) == type(self.molecule_type[_mol].end_type))
                 and (x not in initial_monomers)
                 and (all([len(x) >= len(y) for y in initial_monomers]))
                 and (x[len(self.molecule_type[x].start_type.smiles):] in _mol)]
        if not _mols:
            raise NextIterationException("No hydrolysis")
        product2_smiles = random.choice(_mols)
        product2_mol = self.molecule_type[product2_smiles]

        # Removing molecule
        product1_smiles = _mol[:-(len(product2_smiles) - len(product2_mol.start_type.smiles))] \
                + (TypeA.smiles if isinstance(product2_mol.start_type, TypeB) else TypeB.smiles)

        return {product1_smiles: Molecule(product1_smiles, self.molecule_type[_mol].start_type,
                                          (TypeA if isinstance(product2_mol.start_type, TypeB) else TypeB),
                                          Polymer()),
                product2_smiles: product2_mol}

    def outcome_smiles(self, smiles1: str, smiles2: str) -> Molecule:
        """
        Generates the product of the reaction
        :param smiles1: molecule 1
        :param smiles2: molecule 2
        :return: product
        """

        mol1 = self.molecule_type[smiles1]
        mol2 = self.molecule_type[smiles2]

        if type(mol1.end_type) != type(mol2.start_type):
            # Removes funcional groups
            smiles1 = smiles1
            smiles2 = smiles2
            # Result
            result = smiles1 + smiles2
            start_type = mol1.start_type
            end_type = mol2.end_type
        elif type(mol1.end_type) != type(mol2.end_type):
            # Removes funcional groups
            smiles1 = smiles1
            smiles2 = smiles2
            # Result
            result = smiles1 + smiles2[::-1]
            start_type = mol1.start_type
            end_type = mol2.start_type
        elif type(mol1.start_type) != type(mol2.end_type):
            # Removes funcional groups
            smiles1 = smiles1
            smiles2 = smiles2
            # Result
            result = smiles2 + smiles1
            start_type = mol2.start_type
            end_type = mol1.end_type
        elif type(mol1.start_type) != type(mol2.start_type):
            # Removes funcional groups
            smiles1 = smiles1
            smiles2 = smiles2
            # Result
            result = smiles2[::-1] + smiles1
            start_type = mol2.end_type
            end_type = mol1.end_type

        return Molecule(result, start_type, end_type)

    def reaction(self, smiles1: str, smiles2: str) -> None:
        """
        Reacts 2 molecules
        :param smiles1: smiles of molecule 1
        :param smiles2: smiles of molecule 2
        :return: None
        """

        # Get the molecules
        mol1 = self.molecule_type[smiles1]
        mol2 = self.molecule_type[smiles2]

        # If the molecules do not react with each other
        if mol1.end_type == mol1.start_type == mol2.end_type == mol2.start_type:
            return

        # Check if at least one is a byproduct
        if isinstance(mol1.type, ByProduct) or isinstance(mol2.type, ByProduct):
            product_smiles = self.reverse_reaction(smiles1, smiles2)
            for _smi, _m in product_smiles.items():
                self.update_state(_smi, 1)
                self.molecule_type[_smi] = _m
        # Reacting
        else:
            product_smiles = self.outcome_smiles(smiles1, smiles2)
            # Update molecule_type dictionary
            self.molecule_type[product_smiles.smiles] = product_smiles
            self.molecule_type[ByProduct.smiles] = Molecule(ByProduct.smiles, None, None, ByProduct())
            # Updates the state
            self.update_state(product_smiles.smiles, 1)

        # Updates the state
        self.update_state(smiles1, -1)
        self.update_state(smiles2, -1)

    def reaction_combinations(self) -> list:
        """
        Makes all possible combinations of reactions
        :return: list of possible combinations
        """
        # Takes only the molecules with more than 1 occurency
        molecules = [x for x in self.state.keys() if (self.state[x] > 0)]
        # Makes all possible combinations
        combinations = [x for x in itertools.combinations(molecules, 2) if not (self.molecule_type[x[0]].end_type == self.molecule_type[x[0]].start_type
                                                                                == self.molecule_type[x[1]].end_type == self.molecule_type[x[1]].start_type)]
        return combinations

    def update_reaction_probability(self) -> None:
        """
        Updates the reaction probability dictionary
        :return: None
        """
        # Create pair_dict
        self.pair_dict = {id: pair for id, pair in enumerate(self.reaction_combinations())}

        # Calculate each rate of reaction
        self.rate_of_reaction_dict = {}
        for id, pair in self.pair_dict.items():
            self.rate_of_reaction_dict[id] = self.rate_of_reaction(*pair)

        _sum = sum([x for x in self.rate_of_reaction_dict.values()])

        # Check if _sum is zero
        if _sum == 0:
            raise EndSimulationException("All rates are 0.")

        # Calculate each probability
        self.reaction_probability_dict = {id: value / _sum for id, value in self.rate_of_reaction_dict.items()}

    def choose_reaction(self) -> None:
        """
        Chooses the reaction based on reaction probability
        :return: None
        """

        # Choose random number
        r = random.random()

        # Lists probability of each pair
        keys, rates = zip(*self.reaction_probability_dict.items())

        # Choose a pair
        for u in range(len(keys)):
            sum_u = sum(rates[:u])
            sum_u_plus1 = sum(rates[:u+1])
            type = keys[u]
            if sum_u < r <= sum_u_plus1:
                break

        # Determine the pair
        id = int(keys[u])

        self.chosen_reaction = self.pair_dict[id]
        self.reaction_type = type

    def choose_tau(self) -> None:
        """
        Calculates tau
        :return: None
        """

        # Choose random number
        r = random.random()

        # Calculate tau
        _sum = sum([x for x in self.rate_of_reaction_dict.values()])
        tau = math.log(1/r)/_sum

        self.tau = tau

    def update_state(self, smiles: Union[list, str], increase: int) -> None:
        """
        Removes or adds molecules to state
        :param smiles: molecule to be updated
        :param increase: how much it decreased/increased
        :return: None
        """
        if type(smiles) == list:
            for el in smiles:
                # Adds if it is new
                if el not in self.state.keys():
                    self.state[el] = 1
                else:
                    self.state[el] += increase
                # Deletes from dict if contains 0 molecules
                if self.state[el] == 0:
                    del self.state[el]
        else:
            # Adds if it is new
            if (smiles not in self.state.keys()) and (smiles[::-1] not in self.state.keys()):
                self.state[smiles] = 1
            elif smiles in self.state.keys():
                self.state[smiles] += increase
            else:
                self.state[smiles[::-1]] += increase
            # Deletes from dict if contains 0 molecules
            try:
                if self.state[smiles] == 0:
                    del self.state[smiles]
            except:
                if self.state[smiles[::-1]] == 0:
                    del self.state[smiles[::-1]]

    def simulate_step(self) -> None:
        """
        Simulates the polymerization at a single step. Computes the reaction probability, chooses the reaction,
        calculates tau and react
        :return: None
        """
        self.update_reaction_probability()
        self.choose_reaction()
        self.choose_tau()
        self.reaction(*self.chosen_reaction)

    def simulate_system(self):
        """
        Runs the simulation from t to t_max and plots the results
        :return: None
        """

        # Runs the simulation
        try:

            while self.t <= self.t_max:
                try:
                    self.simulate_step()
                    self.t += self.tau
                    self.time_evolution[self.t] = copy.deepcopy(self.state)
                except NextIterationException as nie:
                    print(nie)

        except EndSimulationException as ese:
            print(ese)

        # Saves the results
        self.save()


    def save(self) -> None:
        """
        Saves the state and the time evolution as a .jsons
        :return: None
        """
        # Saving the state
        state = {}
        state["final_state"] = {k: v for k, v in self.state.items() if v != 0}
        state["input_parameters"] = input_parameters
        with open(output_folder + "Final_State_new.pickle", "wb") as f:
            pickle.dump(state, f, protocol=pickle.HIGHEST_PROTOCOL)
        # Saving the time evolution
        with open(output_folder + 'Time_evolution_new.pickle', 'wb') as f:
            pickle.dump(self.time_evolution, f, protocol=pickle.HIGHEST_PROTOCOL)
        # Saving molecule_type
        with open(output_folder + 'Molecule_type_dict_new.pickle', 'wb') as f:
            pickle.dump(self.molecule_type, f, protocol=pickle.HIGHEST_PROTOCOL)

# Main -----------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":


    # Creating folders
    Paths.create_folder(output_folder)

    # Simulation setup
    initial_system = {mol["molecule"]: number_of_molecules/len(molecules) for mol in molecules}

    # Molecule type
    molecule_type = {}
    for mol in molecules:
        # Smiles of the molecule
        smiles = mol["molecule"]
        # Choosing endgroup type
        if mol["endgroup"] == TypeA.smiles:
            endgroup_type = TypeA()
        else:
            endgroup_type = TypeB()
        # Storing in the dict
        molecule_type[smiles] = Molecule(smiles, endgroup_type, endgroup_type)


    sys = System_FromMonomer(entry_molecules_smiles=initial_system, molecule_type=molecule_type, t=t,
                 t_max=t_max, kinetics=kinetics, V=number_of_molecules/(initial_concentration*Na))

    t0 = datetime.datetime.now()
    # Running simulation
    sys.simulate_system()
    t = datetime.datetime.now()
    print(t-t0)

