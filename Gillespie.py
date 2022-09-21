from operator import eq
import random
import math
import re
from matplotlib import pyplot as plt

class StructuralFormula:
    def __init__(self,chemical,coefficient):
        self.chemical = chemical
        self.coefficient = coefficient
    def getChemical(self):
        return self.chemical
    def getCoefficient(self):
        return self.coefficient

class Reaction:
    def __init__(self,reactants = None,products = None, rate = None,equation = None):
        # Initalization, try to parse a full 'A + B -> C' equation
        if equation is not None:
            self.parseFullRxn(equation)
        else:
            # otherwise is reactant and products are string parse those individually
            if isinstance(reactants, str):
                self.reactants = self.parseRxn(reactants)
            else: 
                #other wise just assume we've been given structural formulas
                self.reactants = reactants

            if isinstance(products, str):
                self.products = self.parseRxn(products)
            else:
                self.products = products

        self.rate = rate

    def parseFullRxn(self,equation):
        # Reactions are of the form 'a+b->c'
        # we first strip to get rid of extra white space
        # than we split along the arrow to get the products and reactants 
        (reactants,products) = equation.strip().split('->')
        self.reactants = self.parseRxn(reactants)
        self.products = self.parseRxn(products)
    def parseRxn(self,eqn):
        # Each side of the equation is a set of possible numbers followed by
        # chemical labels followed sepearted by plus signs
        chemicals = eqn.strip().split('+')
        structuralFormulas = []
        for chemical in chemicals:
            #find the (possible) coeffeient and the label
            # look for any number of digits grouped together
            coeffMatch = re.search('(\d+)',chemical)
            coeff = 1 # default to one
            # if we found a match
            if coeffMatch:
                # set the coefficent equal to it
                coeff = int(coeffMatch.group(1))
            label = ''
            # find the label for the chemical
            labelMatch = re.search('([A-Za-z]+)',chemical)
            if labelMatch:
                label = labelMatch.group(1)
            else:
                # if we can't find one throw and error
                raise Exception("Could not parse chemical equation %s"%chemical)
            
            structuralFormulas.append(StructuralFormula(label,coeff))
        return structuralFormulas

    def getPropensity(self,chemicalList):
        # gets the propensity of a reaction given a list of chemicals
        a = self.rate
        for chemical in self.reactants:
            for i in range(chemical.getCoefficient()):
                a *= chemicalList[chemical.getChemical()]
        return a

    def performReaction(self,chemicalList):
        # repeated code... an opportunity for refactoring
        # for each chemical in the reactants reduce the number by the coefficient
        # and for each chemical in the products increase it by the coefficient 
        for chemical in self.reactants:
            chemicalList[chemical.getChemical()]-= chemical.getCoefficient()
        for chemical in self.products:
            chemicalList[chemical.getChemical()]+= chemical.getCoefficient()
        
    
class Gillespie:
    def __init__(self,initialValues,reactions):
        # initalize the chemical reations with some inital values and a list of 
        # reactions to perform
        self.initialValues = initialValues # save intial values for resetting later
        self.chemicalList = initialValues
        self.reactions = reactions
        self.propensities = [0.0] * len(reactions) # the list of propensities
        self.totalP = 0 # the sum total of propensities
        self.T = 0 # Current time
        self.observer = None # observer to take measurements

    def setObserver(self,obs):
        self.observer = obs 
    
    def reset(self):
        # reset to initial conditions
        self.T = 0
        self.chemicalList = self.initialValues
    
    def getT(self):
        return self.T
    
    def getChemicals(self):
        # if our chemicals are in a list return a tuple of the list
        if isinstance(self.chemicalList, list): 
            return tuple(self.chemicalList)

        # if our chemicals are in a dictionary return a tuple of the values
        if isinstance(self.chemicalList,dict):
            return tuple(self.chemicalList.values())

    def run(self,maxT):
        # run the equation up to a time maxT
        while self.T<maxT:
            self.step()

    def step(self):
        # perform a reaction and update the time
        self.calculatePropensities()
        self.updateTime()
        RxNum = self.chooseReaction()
        self.reactions[RxNum].performReaction(self.chemicalList)
        # if we have an observer take a measurement
        if self.observer is not None:
            self.observer.takeMeasurement(self)
        
    def calculatePropensities(self):
        self.totalP = 0
        # for each reaction ask it for its propensites and then keep a running sum
        for i in range(len(self.reactions)):
            self.propensities[i] = self.reactions[i].getPropensity(self.chemicalList)
            self.totalP += self.propensities[i]
    
    def updateTime(self):
        dt = 1/self.totalP*math.log(1/random.random())
        self.T += dt
    
    def chooseReaction(self):
        # pick a reaction to do weighted by the propensity of that reaction
        r = random.random()
        RxNum = 0
        tally = 0
        for a in self.propensities:
            tally += a
            if tally/self.totalP>r:
                break
            RxNum+=1
        return RxNum

class GillespieObserver:
    # takes measurments and makes plots for a gillespie algorithm
    def __init__(self,fileName):
        # file name will be the name of the data file and the plot
        self.fileName = fileName
        self.file  = open(fileName+".txt",'w')
        self.T = []
        self.chems = []

    def __del__(self):
        self.file.close()

    def takeMeasurement(self,gillespie):
        # write the time to a file and save it in the time list
        self.file.write(str(gillespie.getT())+" ")
        self.T.append(gillespie.getT())
        # write the chemicals to a file and save them to the chemical list
        self.file.write(str(gillespie.getChemicals())+"\n")
        self.chems.append(gillespie.getChemicals())

    def Plot(self):
        # create a time series of the chemical trajectories
        fig,ax = plt.subplots()
        plt.plot(self.T,self.chems)        
        plt.savefig(self.fileName+".png")

if __name__=="__main__":
    # A + B -> C
    # rx0 = Reaction([StructuralFormula(0,1),StructuralFormula(1,1)],[StructuralFormula(2,1)],1)
    # C -> A + B
    # rx1 = Reaction([StructuralFormula(2,1)],[StructuralFormula(0,1),StructuralFormula(1,1)],1)
    # myGillespie = Gillespie([50,50,50],[rx0,rx1])

    rx0 = Reaction(equation='A + B -> C',rate = 1)
    rx1 = Reaction(equation='C->A+B',rate = 1)
    myGillespie = Gillespie({'A':50,'B':50,'C':50},[rx0,rx1])

    myObs = GillespieObserver("data")

    myGillespie.setObserver(myObs)

    myGillespie.run(10)
    myObs.Plot()
    plt.show()
