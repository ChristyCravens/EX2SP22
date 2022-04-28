import numpy as np
import math
from scipy.optimize import fsolve
import random as rnd

#JES broken code in Fluid class
class Fluid():
    def __init__(self, mu=0.00003050735995, rho=997.7733924224):
        '''
        default properties are for water
        :param mu: dynamic viscosity in lb-s/ft^2
        :param rho: density in kg/m^3
        '''
        self.mu= mu # simply make a copy of the value in the argument as a class property
        self.rho= rho # simply make a copy of the value in the argument as a class property
        self.nu= mu/rho # calculate the kinematic viscosity in units of m^2/s

#JES broken code in Node class
class Node():
    def __init__(self, Name='a', Pipes=[], ExtFlow=0):
        '''
        A node in a pipe network.
        :param Name: name of the node
        :param Pipes: a list/array of pipes connected to this node
        :param ExtFlow: any external flow into (+) or out (-) of this node in L/s
        '''
        self.name=Name
        self.pipes=Pipes
        self.extFlow=ExtFlow

    def getNetFlowRate(self):
        '''
        Calculates the net flow rate into this node in L/s
        # :return:
        '''
        Qtot= self.extFlow #$JES MISSING CODE$  #count the external flow first
        for p in self.pipes:
            #retrieves the pipe flow rate (+) if into node (-) if out of node.  see class for pipe.
            Qtot+=p.getFlowIntoNode(self.name)
        return Qtot*0.0353146667

#JES nothing broken in Loop class
class Loop():
    def __init__(self, Name='A', Pipes=[]):
        '''
        Defines a loop in a pipe network.  Note: the pipes must be listed in order.  The traversal of a pipe loop
        will begin at the start node of Pipe[0] and move in the positive direction of that pipe.  Hence, loops
        can be either CW or CCW traversed, depending on which pipe you start with.  Should work fine either way.
        :param Name: name of the loop
        :param Pipes: a list/array of pipes in this loop
        '''
        self.name=Name
        self.pipes=Pipes

    def getLoopHeadLoss(self):
        '''
        Calculates the net head loss as I traverse around the loop, in m of fluid.
        :return:
        '''
        deltaP=0 #initialize to zero
        startNode=self.pipes[0].startNode #begin at the start node of the first pipe
        for p in self.pipes:
            # calculates the head loss in the pipe considering loop traversal and flow directions
            phl=p.getFlowHeadLoss(startNode)
            deltaP+=phl
            startNode=p.endNode if startNode!=p.endNode else p.startNode #move to the next node
        return deltaP

#JES broken code in pipe class
class Pipe():
    def __init__(self, Start='A', End='B', L=100, D=200, r=0.003, fluid=Fluid()):
        '''
        Defines a generic pipe with orientation from lowest letter to highest, alphabetically.
        :param Start: the start node (string)
        :param End: the end node (string)
        :param L: the pipe length in m (float)
        :param D: the pipe diameter in mm (float)
        :param r: the pipe roughness in m  (float)
        :param fluid:  a Fluid object (typically water)
        '''
        # from arguments given in constructor
        self.startNode=min(Start,End) #makes sure to use the lowest letter for startNode
        self.endNode=max(Start,End) #makes sure to use the highest letter for the endNode
        self.length=L
        self.r=r
        self.fluid=fluid #the fluid in the pipe

        # other calculated properties
        self.d=D/1000.0 #diameter in m
        self.relrough = self.r/self.d #calculate relative roughness for easy use later
        self.A=math.pi/4.0*self.d**2 #calculate pipe cross-sectional area for easy use later
        self.Q=10 #working in units of L/s, just an initial guess
        self.vel=self.V()  #calculate the initial velocity of the fluid
        self.reynolds=self.Re() #calculate the initial reynolds number

    def V(self):
        '''
        Calculate average velocity in the pipe for volumetric flow self.Q
        :return:the average velocity in m/s
        '''
        self.vel= ((self.Q/1000) / self.A) #$JES MISSING CODE$  # the average velocity is Q/A (be mindful of units)
        # Q is in L/s, 1000 liters per cubic meter, divide by 1000 gives m^3/s
        return self.vel

    def Re(self):
        '''
        Calculate the reynolds number under current conditions.
        :return:
        '''
        self.reynolds= (self.V()*self.length)/self.fluid.nu #$JES MISSING CODE$ # Re=rho*V*d/nu, be sure to use V() so velocity is updated.
        return self.reynolds

    def FrictionFactor(self):
        """
        This function calculates the friction factor for a pipe based on the
        notion of laminar, turbulent and transitional flow.
        :return: the (Darcy) friction factor
        """
        # update the Reynolds number and make a local variable Re
        Re=self.Re()
        rr=self.relrough
        # to be used for turbulent flow
        def CB():
            # note:  in numpy log is for natural log.  log10 is log base 10.
            cb = lambda f: 1 / (f ** 0.5) + 2.0 * np.log10(rr / 3.7 + 2.51 / (Re * f ** 0.5))
            result = fsolve(cb, (0.01))
            val = cb(result[0])
            return result[0]
        # to be used for laminar flow
        def lam():
            return 64 / Re

        if Re >= 4000:  # true for turbulent flow
            return CB()
        if Re <= 2000:  # true for laminar flow
            return lam()

        # transition flow is ambiguous, so use normal variate weighted by Re
        CBff = CB()
        Lamff = lam()
        # I assume laminar is more accurate when just above 2000 and CB more accurate when just below Re 4000.
        # I will weight the mean appropriately using a linear interpolation.
        mean = Lamff+((Re-2000)/(4000-2000))*(CBff - Lamff)
        sig = 0.2 * mean
        # Now, use normalvariate to put some randomness in the choice
        return rnd.normalvariate(mean, sig)

    def frictionHeadLoss(self):  # calculate headloss through a section of pipe in m of fluid
        '''
        Use the Darcy-Weisbach equation to find the head loss through a section of pipe.
        '''
        g = 9.81  # m/s^2
        ff = self.FrictionFactor()
        hl = ff*((self.length/self.d)) * (self.vel**2/(2*g)) #$JES MISSING CODE$
        return hl

    def getFlowHeadLoss(self, s):
        '''
        Calculate the head loss for the pipe.
        :param s: the node i'm starting with in a traversal of the pipe
        :return: the signed headloss through the pipe in m of fluid
        '''
        #while traversing a loop, if s = startNode I'm traversing in same direction as positive pipe
        nTraverse= 1 if s==self.startNode else -1
        #if flow is positive sense, scalar =1 else =-1
        nFlow=1 if self.Q >= 0 else -1
        return nTraverse*nFlow*self.frictionHeadLoss()

    def Name(self):
        '''
        Gets the pipe name.
        :return:
        '''
        return self.startNode+'-'+self.endNode

    def oContainsNode(self, node):
        #does the pipe connect to the node?
        return self.startNode==node or self.endNode==node

    def printPipeFlowRate(self):
        Q = self.Q*0.0353146667 # Conversion to cfs from L/s
        print('The flow in segment {} is {:0.2f} cfs'.format(self.Name(),Q)) # Convert from SI to English here. Boolean SI=True, set to false for english units, etc.

    def getFlowIntoNode(self, n):
        '''
        determines the flow rate into node n
        :param n: a node object
        :return: +/-Q
        '''
        if n==self.startNode:
            return -self.Q
        return self.Q

#JES broken code in pipe class
class PipeNetwork():
    def __init__(self, Pipes=[], Loops=[], Nodes=[], fluid=Fluid()):
        '''
        The pipe network is built from pipe, node, loop, and fluid objects.
        :param Pipes: a list of pipe objects
        :param Loops: a list of loop objects
        :param Nodes: a list of node objects
        :param fluid: a fluid object
        '''
        self.loops=Loops
        self.nodes=Nodes
        self.Fluid=fluid
        self.pipes=Pipes

    def findFlowRates(self):
        '''
        a method to analyze the pipe network and find the flow rates in each pipe
        given the constraints of: i) no net flow into a node and ii) no net pressure drops in the loops.
        :return: a list of flow rates in the pipes
        '''
        #see how many nodes and loops there are, this is how many equation results I will return
        N=len(self.nodes)+len(self.loops)
        # build an initial guess for flow rates in the pipes.
        # note that I only have 10 pipes, but need 11 variables because of the degenerate node equation at b.
        Q0=np.full(N,14) # CHANGED FROM 10 TO 14; 13 PIPES NOW
        def fn(q):
            """
            This is used as a callback for fsolve.  The mass continuity equations at the nodes and the loop equations
            are functions of the flow rates in the pipes.  Hence, fsolve will search for the roots of these equations
            by varying the flow rates in each pipe.
            :param q: an array of flowrates in the pipes + 1 extra value b/c of node b
            :return: L an array containing flow rates at the nodes and  pressure losses for the loops
            """
            #update the flow rate in each pipe object
            for i in range(len(self.pipes)):
                self.pipes[i].Q= q[i] #$JES MISSING CODE$  # set volumetric flow rate from input argument q
            #calculate the net flow rate for the node objects
            # note:  when flow rates in pipes are correct, the net flow into each node should be zero.
            L= self.getNodeFlowRates() #$JES MISSING CODE$  # call the getNodeFlowRates function of this class
            #calculate the net head loss for the loop objects
            # note: when the flow rates in pipes are correct, the net head loss for each loop should be zero.
            L+= self.getLoopHeadLosses() #$JES MISSING CODE$  # call the getLoopHeadLoss function of this class
            return L
        #using fsolve to find the flow rates
        FR=fsolve(fn,Q0)
        return FR

    def getNodeFlowRates(self):
        #each node object is responsible for calculating its own net flow rate
        qNet=[n.getNetFlowRate() for n in self.nodes]
        return qNet

    def getLoopHeadLosses(self):
        #each loop object is responsible for calculating its own net head loss
        lhl=[l.getLoopHeadLoss() for l in self.loops]
        return lhl

    def getPipe(self, name):
        #returns a pipe object by its name
        for p in self.pipes:
            if name == p.Name():
                return p

    def getNodePipes(self, node):
        #returns a list of pipe objects that are connected to the node object
        l=[]
        for p in self.pipes:
            if p.oContainsNode(node):
                l.append(p)
        return l

    def nodeBuilt(self, node):
        #determines if I have already constructed this node object (by name)
        for n in self.nodes:
            if n.name==node:
                return True
        return False

    def getNode(self, name):
        #returns one of the node objects by name
        for n in self.nodes:
            if n.name==name:
                return n

    def buildNodes(self):
        #automatically create the node objects by looking at the pipe ends
        for p in self.pipes:
            if self.nodeBuilt(p.startNode)==False:
                #instantiate a node object and append it to the list of nodes
                self.nodes.append(Node(p.startNode,self.getNodePipes(p.startNode)))
            if self.nodeBuilt(p.endNode)==False:
                #instantiate a node object and append it to the list of nodes
                self.nodes.append(Node(p.endNode,self.getNodePipes(p.endNode)))

    def printPipeFlowRates(self):
        for p in self.pipes:
            p.printPipeFlowRate()

    def printNetNodeFlows(self):
        for n in self.nodes:
            print('net flow into node {} is {:0.2f}'.format(n.name, n.getNetFlowRate()*0.0353146667)) # convert to English units

    def printLoopHeadLoss(self):
        for l in self.loops:
            print('head loss for loop {} is {:0.2f}'.format(l.name, l.getLoopHeadLoss()))

    def printHeadLoss(self):
        """
        This function finds the head loss in each segment of the pipe system using the frictionHeadLoss function,
        then converts it from meters to inches and prints the head loss in each segment.
        :return:
        """
        # for loop to loop through each segment of pipe and calculate the head loss
        for h in self.pipes:
        # Calculate head loss in inches of water, converting from meters to inches
            hL = h.frictionHeadLoss()*(1/0.0254) # m*(1/0.0254)=inches
            # print each segment formatting the name of the segment and the head loss from the above equation.
            print('The head loss in segment {} is {:0.2f} inches'.format(h.Name(),hL))

def Pressure():
    """
    This function finds the pressures in each segment of pipe and returns it as a list called "Pressures"
    :return: list of pressures
    """
    # Initial values to use in loop for all nodes of the pipe system
    # Specific weight
    Y = 62.3 * (1/1728) # Converting from lb/ft^3 to lb/in^3
    y = [19200, 19200, 19200, 19200, 9600, 9600, 9600, 0, 0, 0] # Distance from node h, in inches
    k = 0
    # Initial pressure
    P0=80 # psi at node h
    Pressures=[]
    # for loop that goes through the given values above and stores each pressure value accordingly
    for x in range(10):
        # This equation solves for pressure in each pipe, using the head loss from above
        P = Y * (y[x]-k)+P0
        Pressures.append(P)
        x += 1
    return Pressures

def main():
    '''
    This program analyzes flows in a given pipe network based on the following:
    1. The pipe segments are named by their endpoint node names:  e.g., a-b, b-e, etc.
    2. Flow from the lower letter to the higher letter of a pipe is considered positive.
    3. Pressure decreases in the direction of flow through a pipe.
    4. At each node in the pipe network, mass is conserved.
    5. For any loop in the pipe network, the pressure loss is zero
    Approach to analyzing the pipe network:
    Step 1: build a pipe network object that contains pipe, node, loop and fluid objects
    Step 2: calculate the flow rates in each pipe using fsolve
    Step 3: output results
    Step 4: check results against expected properties of zero head loss around a loop and mass conservation at nodes.
    :return:
    '''
    #instantiate a Fluid object to define the working fluid as water
    water= Fluid()
    roughnessA = 0.00025908  # in m, 12" & 16" pipes
    roughnessB = 0.0009144 # in m, 18" & 24" pipes

    #instantiate a new PipeNetwork object
    Pipes=[]
    Loops=[]
    PN= PipeNetwork(Pipes, Loops)#$JES MISSING CODE$  #
    P= Pipe()
    #add Pipe objects to the pipe network (see constructor for Pipe class)
    PN.pipes.append(Pipe('a','b',304.8, 457.2, roughnessB, water))
    PN.pipes.append(Pipe('a','h',487.68, 609.6, roughnessB, water))
    PN.pipes.append(Pipe('b','e',243.84, 406.4, roughnessA, water))
    PN.pipes.append(Pipe('b','c',152.4, 457.2, roughnessB, water))
    PN.pipes.append(Pipe('c','d',152.4, 457.2, roughnessB, water))
    PN.pipes.append(Pipe('c','f',243.84, 406.4, roughnessA, water))
    PN.pipes.append(Pipe('d','g',243.84, 406.4, roughnessA, water))
    PN.pipes.append(Pipe('e','f',152.4, 304.8, roughnessA, water))
    PN.pipes.append(Pipe('e','i',243.84, 457.2, roughnessB, water))
    PN.pipes.append(Pipe('f','g',152.4, 304.8, roughnessA, water))
    PN.pipes.append(Pipe('g','j',243.84, 457.2, roughnessB, water))
    PN.pipes.append(Pipe('h','i',304.8, 609.6, roughnessB, water))
    PN.pipes.append(Pipe('i','j',304.8, 609.6, roughnessB, water))

    #add Node objects to the pipe network by calling buildNodes method of PN object
    PN.buildNodes()

    #update the external flow of certain nodes
    PN.getNode('d').extFlow=-2*28.316847 # conversion from cfs to L/s
    PN.getNode('e').extFlow=-3*28.316847
    PN.getNode('f').extFlow=-5*28.316847
    PN.getNode('h').extFlow=10*28.316847

    #add Loop objects to the pipe network
    PN.loops.append(Loop('A',[PN.getPipe('a-b'), PN.getPipe('b-e'),PN.getPipe('e-i'), PN.getPipe('h-i'), PN.getPipe('a-h')]))
    PN.loops.append(Loop('B',[PN.getPipe('b-c'), PN.getPipe('c-f'),PN.getPipe('e-f'), PN.getPipe('b-e')]))
    PN.loops.append(Loop('C',[PN.getPipe('c-d'), PN.getPipe('d-g'),PN.getPipe('f-g'), PN.getPipe('c-f')]))
    PN.loops.append(Loop('D', [PN.getPipe('e-f'), PN.getPipe('f-g'), PN.getPipe('g-j'), PN.getPipe('i-j'), PN.getPipe('e-i')]))

    #call the findFlowRates method of the PN (a PipeNetwork object)
    PN.findFlowRates()

    #get output
    print()
    print('Flow in each segment of pipe:')
    PN.printPipeFlowRates()
    print()
    print('Check node flows:')
    PN.printNetNodeFlows()
    print()
    print('Check loop head loss:')
    PN.printLoopHeadLoss()
    print()
    print('Head loss for segments:')
    PN.printHeadLoss()
    print()
    print('Pressure in each node:')
    P = Pressure()
    nodes = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
    for x in range(10):
        print('The pressure in node {} is {:0.2f} psi'.format(nodes[x], P[x]))

main()