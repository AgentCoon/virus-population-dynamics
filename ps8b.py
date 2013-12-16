# Problem Set: Simulating the Spread of Disease and Virus Population Dynamics 

import numpy
import random
import pylab

''' 
Begin helper code
'''

class NoChildException(Exception):
    """
    NoChildException is raised by the reproduce() method in the SimpleVirus
    and ResistantVirus classes to indicate that a virus particle does not
    reproduce. You can use NoChildException as is, you do not need to
    modify/add any code.
    """

'''
End helper code
'''

#
# PROBLEM 1
#
class SimpleVirus(object):

    """
    Representation of a simple virus (does not model drug effects/resistance).
    """
    def __init__(self, maxBirthProb, clearProb):
        """
        Initialize a SimpleVirus instance, saves all parameters as attributes
        of the instance.        
        maxBirthProb: Maximum reproduction probability (a float between 0-1)        
        clearProb: Maximum clearance probability (a float between 0-1).
        """
        self.maxBirthProb = maxBirthProb
        self.clearProb = clearProb


    def getMaxBirthProb(self):
        """
        Returns the max birth probability.
        """
        return self.maxBirthProb

    def getClearProb(self):
        """
        Returns the clear probability.
        """
        return self.clearProb

    def doesClear(self):
        """ Stochastically determines whether this virus particle is cleared from the
        patient's body at a time step. 
        returns: True with probability self.getClearProb and otherwise returns
        False.
        """
        clear = random.random()
        return clear < self.getClearProb()

    def reproduce(self, popDensity):
        """
        Stochastically determines whether this virus particle reproduces at a
        time step. Called by the update() method in the Patient and
        TreatedPatient classes. The virus particle reproduces with probability
        self.getMaxBirthProb * (1 - popDensity).
        
        If this virus particle reproduces, then reproduce() creates and returns
        the instance of the offspring SimpleVirus (which has the same
        maxBirthProb and clearProb values as its parent).         

        popDensity: the population density (a float), defined as the current
        virus population divided by the maximum population.         
        
        returns: a new instance of the SimpleVirus class representing the
        offspring of this virus particle. The child should have the same
        maxBirthProb and clearProb values as this virus. Raises a
        NoChildException if this virus particle does not reproduce.               
        """

        ran = random.random()
        if ran < self.getMaxBirthProb() * (1 - popDensity):
            child = SimpleVirus(self.getMaxBirthProb(), self.getClearProb())
            return child
        else:
            raise NoChildException()


class Patient(object):
    """
    Representation of a simplified patient. The patient does not take any drugs
    and his/her virus populations have no drug resistance.
    """    
    def __init__(self, viruses, maxPop):
        """
        Initialization function, saves the viruses and maxPop parameters as
        attributes.

        viruses: the list representing the virus population (a list of
        SimpleVirus instances)

        maxPop: the maximum virus population for this patient (an integer)
        """
        self.viruses = viruses
        self.maxPop = maxPop
        
    def getViruses(self):
        """
        Returns the viruses in this Patient.
        """
        return self.viruses

    def getMaxPop(self):
        """
        Returns the max population.
        """
        return self.maxPop

    def getTotalPop(self):
        """
        Gets the size of the current total virus population. 
        returns: The total virus population (an integer)
        """
        return len(self.viruses)   

    def update(self):
        """
        Update the state of the virus population in this patient for a single
        time step. update() should execute the following steps in this order:
        
        - Determine whether each virus particle survives and updates the list
        of virus particles accordingly.   
        
        - The current population density is calculated. This population density
          value is used until the next call to update() 
        
        - Based on this value of population density, determine whether each 
          virus particle should reproduce and add offspring virus particles to 
          the list of viruses in this patient.                    

        returns: The total virus population at the end of the update (an
        integer)
        """
        virusList = []
        for v in self.viruses:
            if v.doesClear() == False:
                virusList.append(v)
        self.viruses = virusList
        popDen = (1.0 * self.getTotalPop()) / self.getMaxPop()

        for v in self.viruses[:]:
            try:
                child = v.reproduce(popDen)
                self.viruses.append(child)
            except NoChildException:
                pass            
        return self.getTotalPop()
        
#
# PROBLEM 2
#
def simulationWithoutDrug(numViruses, maxPop, maxBirthProb, clearProb,
                          numTrials):
    """
    Run the simulation and plot the graph for problem 2 (no drugs are used,
    viruses do not have any drug resistance).    
    For each of numTrials trial, instantiates a patient, runs a simulation
    for 300 timesteps, and plots the average virus population size as a
    function of time.

    numViruses: number of SimpleVirus to create for patient (an integer)
    maxPop: maximum virus population for patient (an integer)
    maxBirthProb: Maximum reproduction probability (a float between 0-1)        
    clearProb: Maximum clearance probability (a float between 0-1)
    numTrials: number of simulation runs to execute (an integer)
    """
    virusList = []
    progressList = []
    for x in range(0, 300):
        progressList.append(0.0)  

    for v in range(0, numViruses):
        flu = SimpleVirus(maxBirthProb, clearProb)
        virusList.append(flu)  

    for t in range(0, numTrials):    
        jacek = Patient(virusList, maxPop)
        for a in range(0, 300):
            jacek.update()
            progressList[a] += jacek.getTotalPop()

    for a in range(0, 300):
        progressList[a] /= float(numTrials)

    pylab.plot(range(0, 300), progressList, label = "Virus population")
    pylab.title("Mean number of viruses over time")
    pylab.xlabel("Time steps")
    pylab.ylabel("Number of virus particles")
    pylab.legend(loc = "best")
    pylab.show()

# Uncomment to see this function in action:
#simulationWithoutDrug(15, 100, 0.5, 0.2, 100)

#
# PROBLEM 3
#
class ResistantVirus(SimpleVirus):
    """
    Representation of a virus which can have drug resistance.
    """   
    def __init__(self, maxBirthProb, clearProb, resistances, mutProb):
        """
        Initialize a ResistantVirus instance, saves all parameters as attributes
        of the instance.

        maxBirthProb: Maximum reproduction probability (a float between 0-1)       

        clearProb: Maximum clearance probability (a float between 0-1).

        resistances: A dictionary of drug names (strings) mapping to the state
        of this virus particle's resistance (either True or False) to each drug.
        e.g. {'guttagonol':False, 'srinol':False}, means that this virus
        particle is resistant to neither guttagonol nor srinol.

        mutProb: Mutation probability for this virus particle (a float). This is
        the probability of the offspring acquiring or losing resistance to a drug.
        """
        SimpleVirus.__init__(self, maxBirthProb, clearProb)
        self.resistances = resistances
        self.mutProb = mutProb

    def getResistances(self):
        """
        Returns the resistances for this virus.
        """
        return self.resistances

    def getMutProb(self):
        """
        Returns the mutation probability for this virus.
        """
        return self.mutProb

    def isResistantTo(self, drug):
        """
        Get the state of this virus particle's resistance to a drug. This method
        is called by getResistPop() in TreatedPatient to determine how many virus
        particles have resistance to a drug.       

        drug: The drug (a string)

        returns: True if this virus instance is resistant to the drug, False
        otherwise.
        """
        return self.resistances.get(drug, False)

    def reproduce(self, popDensity, activeDrugs):
        """
        Stochastically determines whether this virus particle reproduces at a
        time step. Called by the update() method in the TreatedPatient class.

        A virus particle will only reproduce if it is resistant to ALL the drugs
        in the activeDrugs list. For example, if there are 2 drugs in the
        activeDrugs list, and the virus particle is resistant to 1 or no drugs,
        then it will NOT reproduce.

        Hence, if the virus is resistant to all drugs
        in activeDrugs, then the virus reproduces with probability:      

        self.getMaxBirthProb * (1 - popDensity).                       

        If this virus particle reproduces, then reproduce() creates and returns
        the instance of the offspring ResistantVirus (which has the same
        maxBirthProb, clearProb, and mutProb values as its parent). The offspring virus
        will have the same maxBirthProb, clearProb, and mutProb as the parent.

        For each drug resistance trait of the virus (i.e. each key of
        self.resistances), the offspring has probability 1-mutProb of
        inheriting that resistance trait from the parent, and probability
        mutProb of switching that resistance trait in the offspring.       

        For example, if a virus particle is resistant to guttagonol but not
        srinol, and self.mutProb is 0.1, then there is a 10% chance that
        that the offspring will lose resistance to guttagonol and a 90%
        chance that the offspring will be resistant to guttagonol.
        There is also a 10% chance that the offspring will gain resistance to
        srinol and a 90% chance that the offspring will not be resistant to
        srinol.

        popDensity: the population density (a float), defined as the current
        virus population divided by the maximum population       

        activeDrugs: a list of the drug names acting on this virus particle
        (a list of strings).

        returns: a new instance of the ResistantVirus class representing the
        offspring of this virus particle. The child should have the same
        maxBirthProb and clearProb values as this virus. Raises a
        NoChildException if this virus particle does not reproduce.
        """
        res = {}
        ran = random.random()
        for d in activeDrugs:
            if self.isResistantTo(d) == False:
                raise NoChildException()

        for k in self.resistances:   
            if self.resistances[k] == True:
                if random.random() <= 1 - self.getMutProb():
                    res[k] = True
                else:
                    res[k] = False
            elif self.resistances[k] == False:
                if random.random() <= 1 - self.getMutProb():
                    res[k] = False
                else:
                    res[k] = True
                
        if ran < self.getMaxBirthProb() * (1 - popDensity):
            return ResistantVirus(self.getMaxBirthProb(), self.getClearProb(), res, self.getMutProb())
        else:
            raise NoChildException()


class TreatedPatient(Patient):
    """
    Representation of a patient. The patient is able to take drugs and his/her
    virus population can acquire resistance to the drugs he/she takes.
    """
    def __init__(self, viruses, maxPop):
        """
        Initialization function, saves the viruses and maxPop parameters as
        attributes. Also initializes the list of drugs being administered
        (which should initially include no drugs).              

        viruses: The list representing the virus population (a list of
        virus instances)

        maxPop: The  maximum virus population for this patient (an integer)
        """
        Patient.__init__(self, viruses, maxPop)
        drugs = []
        self.drugs = drugs
        
    def addPrescription(self, newDrug):
        """
        Administer a drug to this patient. After a prescription is added, the
        drug acts on the virus population for all subsequent time steps. If the
        newDrug is already prescribed to this patient, the method has no effect.

        newDrug: The name of the drug to administer to the patient (a string).

        postcondition: The list of drugs being administered to a patient is updated
        """
        if newDrug not in self.drugs:
            self.drugs.append(newDrug)

    def getPrescriptions(self):
        """
        Returns the drugs that are being administered to this patient.

        returns: The list of drug names (strings) being administered to this
        patient.
        """
        return self.drugs

    def getResistPop(self, drugResist):
        """
        Get the population of virus particles resistant to the drugs listed in
        drugResist.       

        drugResist: Which drug resistances to include in the population (a list
        of strings - e.g. ['guttagonol'] or ['guttagonol', 'srinol'])

        returns: The population of viruses (an integer) with resistances to all
        drugs in the drugResist list.
        """
        virusNr = 0
        viruses = self.getViruses()
        
        for v in viruses:
            count = 0
            for d in drugResist:
                if v.isResistantTo(d) == True:
                    count += 1
                else:
                    break
            if count == len(drugResist):
                virusNr += 1
        return virusNr
                
    def update(self):
        """
        Update the state of the virus population in this patient for a single
        time step. update() should execute these actions in order:

        - Determine whether each virus particle survives and update the list of
          virus particles accordingly

        - The current population density is calculated. This population density
          value is used until the next call to update().

        - Based on this value of population density, determine whether each 
          virus particle should reproduce and add offspring virus particles to 
          the list of viruses in this patient.
          The list of drugs being administered should be accounted for in the
          determination of whether each virus particle reproduces.

        returns: The total virus population at the end of the update (an
        integer)
        """
        virusList = []
       
        for v in self.viruses:
            if v.doesClear() == False:
                virusList.append(v)
        self.viruses = virusList
        popDen = (1.0 * self.getTotalPop()) / self.getMaxPop()

        for v in self.viruses[:]:
            try:
                child = v.reproduce(popDen, self.getPrescriptions())
                self.viruses.append(child)
            except NoChildException:
                pass        
        return self.getTotalPop()

#
# PROBLEM 4
#
def simulationWithDrug(numViruses, maxPop, maxBirthProb, clearProb, resistances,
                       mutProb, numTrials):
    """
    Runs simulations and plots graphs for problem 4.

    For each of numTrials trials, instantiates a patient, runs a simulation for
    150 timesteps, adds guttagonol, and runs the simulation for an additional
    150 timesteps.  At the end plots the average virus population size
    (for both the total virus population and the guttagonol-resistant virus
    population) as a function of time.

    numViruses: number of ResistantVirus to create for patient (an integer)
    maxPop: maximum virus population for patient (an integer)
    maxBirthProb: Maximum reproduction probability (a float between 0-1)        
    clearProb: maximum clearance probability (a float between 0-1)
    resistances: a dictionary of drugs that each ResistantVirus is resistant to
                 (e.g., {'guttagonol': False})
    mutProb: mutation probability for each ResistantVirus particle
             (a float between 0-1). 
    numTrials: number of simulation runs to execute (an integer)
    
    """
    virusList = []
    
    flu = ResistantVirus(maxBirthProb, clearProb, resistances, mutProb)
    for i in range(0, numViruses):
        virusList.append(flu)
    
    progressList = []
    progressList2 = []
    for x in range(0, 300):
        progressList.append(0.0)
        progressList2.append(0.0)
    
    cured = 0
    for i in range(0, numTrials):
        jacek = TreatedPatient(virusList, maxPop)
        for a in range(0, 300):
            if a == 150:
                jacek.addPrescription("guttagonol")
            jacek.update()
            progressList[a] += jacek.getTotalPop()
            progressList2[a] += jacek.getResistPop(["guttagonol"])
        if jacek.getTotalPop() <= 50:
            cured += 1

    for i in range(0, 300):
        progressList[i] /= float(numTrials)
        progressList2[i] /= float(numTrials)
        
    pylab.plot(range(0, 300), progressList, label="Total virus population")
    pylab.plot(range(0, 300), progressList2, label="Resistant virus population")
    pylab.title("Influence of a drug treatment on a virus population")
    pylab.xlabel("Time steps")
    pylab.ylabel("Number of virus particles")
    pylab.legend(loc="best")
    pylab.grid()
    pylab.show()

# Uncomment to see this function in action:
#example1:
#simulationWithDrug(10, 150, 0.1, 0.05, {"guttagonol": True}, 0.1, 50)
#example2:
#simulationWithDrug(15, 200, 0.06, 0.05, {"guttagonol": False, "grimpex" : True}, 0.05, 50)


def simulationWithDrugs(numViruses, maxPop, maxBirthProb, clearProb, resistances,
                       mutProb, numTrials):
    """
    For each of numTrials trials, instantiates a patient, runs a simulation for
    150 timesteps, adds guttagonol, runs the simulation for an additional
    75 timesteps, adds grimpex, and runs the simulation for 150 time steps.
    At the end plots the average virus population size as a function of time.

    numViruses: number of ResistantVirus to create for patient (an integer)
    maxPop: maximum virus population for patient (an integer)
    maxBirthProb: Maximum reproduction probability (a float between 0-1)        
    clearProb: maximum clearance probability (a float between 0-1)
    resistances: a dictionary of drugs that each ResistantVirus is resistant to
                 (e.g., {'guttagonol': False})
    mutProb: mutation probability for each ResistantVirus particle
             (a float between 0-1). 
    numTrials: number of simulation runs to execute (an integer)
    """
    guttagonolStep = 150
    grimpexStep = 225
    virusList = []
    progressList = []
    virus = ResistantVirus(maxBirthProb, clearProb, resistances, mutProb)

    for i in range(0, numViruses):
        virusList.append(virus)
    for i in range(0, 375):
        progressList.append(0.0)
    
    cured = 0    
    for i in range(0, numTrials):
        jacek = TreatedPatient(virusList, maxPop)
        for a in range(0, 375):
            if a == guttagonolStep:
                jacek.addPrescription("guttagonol")
            if a == grimpexStep:
                jacek.addPrescription("grimpex")
            jacek.update()
            progressList[a] += jacek.getTotalPop()
        if jacek.getTotalPop() <= 50:
            cured += 1

    for i in range(0, 375):
        progressList[i] /= float(numTrials)
        
    pylab.plot(range(0, 375), progressList, label="Total virus population")
    pylab.title("Influence of a drug treatment on a virus population")
    pylab.xlabel("Time steps")
    pylab.ylabel("Number of virus particles")
    pylab.legend(loc="best")
    pylab.grid()
    pylab.show()

# Uncomment below to see this function at work:
#simulationWithDrugs(15, 300, 0.2, 0.15, {"guttagonol" : False, "grimpex" : False}, 0.01, 100)
